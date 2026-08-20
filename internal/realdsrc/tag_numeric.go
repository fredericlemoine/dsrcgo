// Numeric tag field encoding, ported from TagTokenizerEncoder/Decoder's
// numeric-field methods (src/TagModeler.cpp:StoreFields/EncodeNextFields/
// StoreNumericField and ReadFields/DecodeNextFields/ReadNumericField).
//
// Implemented: DeltaConst (a plain incrementing-by-a-fixed-amount counter —
// zero bits per record beyond the first), and DeltaVar/ValueVar, each
// either fixed-bit-width or (when the value/delta range is small enough
// that real dsrc's own var_stat_encode heuristic picks it) Huffman-coded
// via the same real-format tree machinery DNA and quality use.
//
// Not implemented: ValueRle/DeltaRle (run-length variants, selected when a
// field's value or delta repeats in long runs — see tag_analyze.go's
// tryRleVal/tryRleDelta, which are still computed correctly so a field
// that needs one of these schemes is detected and reported, not silently
// mis-encoded).
package realdsrc

import (
	"fmt"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/huffman"
)

const huffmanGlobalCapacity = 512 // Field::HUF_GLOBAL_SIZE

// numericCoder holds per-field encode/decode-time state for a numeric
// field: the running "previous value" cursor plus, when var_stat_encode
// applies, the Huffman codes (encode side) or tree (decode side).
type numericCoder struct {
	prevValue int32

	huffmanCodes  []huffman.Code   // encode side, indexed by (value - min)
	huffmanDecode *realHuffmanTree // decode side
}

// storeNumericFieldHeader mirrors StoreFields's numeric-field branch.
func storeNumericFieldHeader(w *bitio.Writer, f *tagField) (*numericCoder, error) {
	if f.numericScheme == numSchemeValueRle || f.numericScheme == numSchemeDeltaRle {
		return nil, fmt.Errorf("realdsrc: tag field needs an RLE numeric scheme (id %d), not yet implemented", f.numericScheme)
	}

	w.PutByte(1) // is_numeric
	w.PutByte(f.numericScheme)
	w.PutWord(uint32(f.minValue))
	w.PutWord(uint32(f.maxValue))

	nc := &numericCoder{}

	switch f.numericScheme {
	case numSchemeDeltaConst:
		w.PutWord(uint32(f.minDelta))
		w.PutWord(uint32(f.maxDelta))
	case numSchemeDeltaVar:
		w.PutWord(uint32(f.minDelta))
		w.PutWord(uint32(f.maxDelta))
		w.PutByte(boolBit8(f.varStatEncode))
		if f.varStatEncode {
			tree, codes := buildNumericHuffman(f.deltaValues, f.minDelta, f.maxDelta)
			storeRealHuffmanTree(w, tree)
			nc.huffmanCodes = codes
		}
	case numSchemeValueVar:
		w.PutByte(boolBit8(f.varStatEncode))
		if f.varStatEncode {
			tree, codes := buildNumericHuffman(f.numValues, f.minValue, f.maxValue)
			storeRealHuffmanTree(w, tree)
			nc.huffmanCodes = codes
		}
	}
	return nc, nil
}

func buildNumericHuffman(freq map[int32]uint32, min, max int32) (*huffman.Tree, []huffman.Code) {
	diff := int(max-min) + 1
	t := huffman.NewTree()
	t.Restart(diff)
	for j := 0; j < diff; j++ {
		t.Insert(freq[min+int32(j)])
	}
	codes := t.Complete()
	return t, codes
}

func boolBit8(b bool) byte {
	if b {
		return 1
	}
	return 0
}

// encodeNumericField mirrors TagTokenizerEncoder::StoreNumericField.
// recordIdx is the encoder's own 0-based record counter (real records
// only, no redundant first pass — unlike the analyzer's).
func encodeNumericField(w *bitio.Writer, f *tagField, nc *numericCoder, recordIdx int, curValue int32) {
	if recordIdx == 0 {
		w.PutBits(uint32(curValue-f.minValue), f.noOfBitsPerValue)
		nc.prevValue = curValue
		return
	}

	switch f.numericScheme {
	case numSchemeDeltaConst:
		// nothing stored: value == prev + minDelta, always.
	case numSchemeDeltaVar:
		toStore := curValue - nc.prevValue - f.minDelta
		if f.varStatEncode {
			c := nc.huffmanCodes[toStore]
			w.PutBits(c.Code, uint(c.Len))
		} else {
			w.PutBits(uint32(toStore), f.noOfBitsPerNum)
		}
	case numSchemeValueVar:
		toStore := curValue - f.minValue
		if f.varStatEncode {
			c := nc.huffmanCodes[toStore]
			w.PutBits(c.Code, uint(c.Len))
		} else {
			w.PutBits(uint32(toStore), f.noOfBitsPerNum)
		}
	}
	nc.prevValue = curValue
}

// loadNumericFieldHeader mirrors ReadFields's numeric-field branch.
func loadNumericFieldHeader(r *bitio.Reader) (*tagField, *numericCoder, error) {
	f := &tagField{isNumeric: true}
	f.numericScheme = r.GetByte()
	f.minValue = int32(r.GetWord())
	f.maxValue = int32(r.GetWord())
	f.noOfBitsPerValue = bitLength(uint32(f.maxValue - f.minValue))

	nc := &numericCoder{}

	switch f.numericScheme {
	case numSchemeDeltaConst, numSchemeDeltaVar:
		f.minDelta = int32(r.GetWord())
		f.maxDelta = int32(r.GetWord())
		f.noOfBitsPerNum = bitLength(uint32(f.maxDelta - f.minDelta))
		f.isDeltaCoding = true
		f.isDeltaConst = f.numericScheme == numSchemeDeltaConst
		if f.numericScheme == numSchemeDeltaVar {
			f.varStatEncode = r.GetByte() != 0
			if f.varStatEncode {
				nc.huffmanDecode = loadRealHuffmanTree(r)
			}
		}
	case numSchemeValueVar:
		f.noOfBitsPerNum = f.noOfBitsPerValue
		f.varStatEncode = r.GetByte() != 0
		if f.varStatEncode {
			nc.huffmanDecode = loadRealHuffmanTree(r)
		}
	case numSchemeValueRle, numSchemeDeltaRle:
		return nil, nil, fmt.Errorf("realdsrc: tag field uses an RLE numeric scheme (id %d), not yet implemented", f.numericScheme)
	default:
		return nil, nil, fmt.Errorf("realdsrc: unrecognized tag numeric scheme %d", f.numericScheme)
	}
	return f, nc, nil
}

// decodeNumericField mirrors TagTokenizerDecoder::ReadNumericField.
func decodeNumericField(r *bitio.Reader, f *tagField, nc *numericCoder, recordIdx int) int32 {
	if recordIdx == 0 {
		v := int32(r.GetBits(f.noOfBitsPerValue)) + f.minValue
		nc.prevValue = v
		return v
	}

	var v int32
	switch f.numericScheme {
	case numSchemeDeltaConst:
		v = nc.prevValue + f.minDelta
	case numSchemeDeltaVar, numSchemeValueVar:
		var raw int32
		if f.varStatEncode {
			raw = nc.huffmanDecode.DecodeSymbol(r)
		} else {
			raw = int32(r.GetBits(f.noOfBitsPerNum))
		}
		if f.numericScheme == numSchemeDeltaVar {
			v = raw + nc.prevValue + f.minDelta
		} else {
			v = raw + f.minValue
		}
	}
	nc.prevValue = v
	return v
}
