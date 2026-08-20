// Text tag field encoding, ported from TagTokenizerEncoder/Decoder's
// text-field methods (src/TagModeler.cpp:StoreFields/EncodeNextFields and
// ReadFields/DecodeNextFields).
//
// A text field's first-record bytes double as a per-position "hamming
// mask": position k is transmitted only once (in the header) if every
// record shares the same byte there, otherwise a real-format Huffman tree
// is built for that position from all 256 possible byte values (most with
// zero frequency, naturally pruned — see huffman.Tree.Complete). Unlike
// DNA/quality's Huffman paths, no separate presence bitmap or dense
// remapping is needed: all 256 byte values are always inserted (many at
// frequency 0), so a tree's leaf id already equals the raw byte value.
// Positions beyond maxFieldStatLen (128) share one catch-all tree,
// indexed at position maxFieldStatLen; trees/codes/decodeTrees below are
// indexed directly by clamped position (with a nil/absent entry at
// constant positions), matching upstream's own direct indexing rather
// than compacting away the unused slots.
package realdsrc

import (
	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/huffman"
)

type textCoder struct {
	trees []*huffman.Tree  // encode side, index = min(position, maxFieldStatLen); nil where constant
	codes [][]huffman.Code // same indexing

	decodeTrees []*realHuffmanTree // decode side, same indexing
}

// storeTextFieldHeader mirrors StoreFields's text-field branch.
func storeTextFieldHeader(w *bitio.Writer, f *tagField) *textCoder {
	w.PutByte(boolBit8(f.isLenConstant))
	w.PutWord(uint32(f.len))
	w.PutWord(uint32(f.maxLen))
	w.PutWord(uint32(f.minLen))
	for _, b := range f.data {
		w.PutByte(b)
	}
	for j := 0; j < f.len; j++ {
		w.PutBit(boolBit(f.hamMask[j]))
	}
	w.AlignByte()

	slots := maxFieldStatLen + 1
	tc := &textCoder{trees: make([]*huffman.Tree, slots), codes: make([][]huffman.Code, slots)}

	buildAndStore := func(idx int) {
		tr := huffman.NewTree()
		tr.Restart(256)
		for k := 0; k < 256; k++ {
			var freq uint32
			if idx < len(f.charFreq) {
				freq = f.charFreq[idx][k]
			}
			tr.Insert(freq)
		}
		tc.codes[idx] = tr.Complete()
		tc.trees[idx] = tr
		storeRealHuffmanTree(w, tr)
	}

	limit := f.maxLen
	if limit > maxFieldStatLen {
		limit = maxFieldStatLen
	}
	for j := 0; j < limit; j++ {
		if j >= f.len || !f.hamMask[j] {
			buildAndStore(j)
		}
	}
	if f.maxLen >= maxFieldStatLen {
		buildAndStore(maxFieldStatLen)
	}
	return tc
}

// encodeTextField mirrors TagTokenizerEncoder::EncodeNextFields's
// text-field branch.
func encodeTextField(w *bitio.Writer, f *tagField, tc *textCoder, value []byte) {
	if !f.isLenConstant {
		w.PutBits(uint32(len(value)-f.minLen), f.noOfBitsPerLen)
	}
	for j, b := range value {
		if j < f.len && f.hamMask[j] {
			continue
		}
		idx := j
		if idx > maxFieldStatLen {
			idx = maxFieldStatLen
		}
		c := tc.codes[idx][b]
		w.PutBits(c.Code, uint(c.Len))
	}
}

// loadTextFieldHeader mirrors ReadFields's text-field branch.
func loadTextFieldHeader(r *bitio.Reader) (*tagField, *textCoder) {
	f := &tagField{}
	f.isLenConstant = r.GetByte() != 0
	f.len = int(r.GetWord())
	f.maxLen = int(r.GetWord())
	f.minLen = int(r.GetWord())
	f.noOfBitsPerLen = bitLength(uint32(f.maxLen - f.minLen))

	f.data = make([]byte, f.len)
	for i := range f.data {
		f.data[i] = r.GetByte()
	}
	f.hamMask = make([]bool, f.len)
	for j := range f.hamMask {
		f.hamMask[j] = r.GetBit() != 0
	}
	r.AlignByte()

	slots := maxFieldStatLen + 1
	tc := &textCoder{decodeTrees: make([]*realHuffmanTree, slots)}

	limit := f.maxLen
	if limit > maxFieldStatLen {
		limit = maxFieldStatLen
	}
	for j := 0; j < limit; j++ {
		if j >= f.len || !f.hamMask[j] {
			tc.decodeTrees[j] = loadRealHuffmanTree(r)
		}
	}
	if f.maxLen >= maxFieldStatLen {
		tc.decodeTrees[maxFieldStatLen] = loadRealHuffmanTree(r)
	}
	return f, tc
}

// decodeTextField mirrors TagTokenizerDecoder::DecodeNextFields's
// text-field branch.
func decodeTextField(r *bitio.Reader, f *tagField, tc *textCoder) []byte {
	fieldLen := f.len
	if !f.isLenConstant {
		fieldLen = int(r.GetBits(f.noOfBitsPerLen)) + f.minLen
	}

	out := make([]byte, fieldLen)
	for k := 0; k < fieldLen; k++ {
		if k < f.len && f.hamMask[k] {
			out[k] = f.data[k]
			continue
		}
		idx := k
		if idx > maxFieldStatLen {
			idx = maxFieldStatLen
		}
		out[k] = byte(tc.decodeTrees[idx].DecodeSymbol(r))
	}
	return out
}
