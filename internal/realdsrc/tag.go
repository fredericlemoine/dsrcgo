// Tag (FASTQ header/ID) encoding entry points, ported from
// TagTokenizerEncoder/TagTokenizerDecoder's top-level methods
// (src/TagModeler.cpp:StartEncoding/StoreFields/EncodeNextFields/
// FinishEncoding and StartDecoding/ReadFields/DecodeNextFields/
// FinishDecoding).
//
// Not implemented: mixed formatting (tags that don't all tokenize to the
// same field count and separators fall back to TagRawEncoder/
// TagRawDecoder upstream — a whole-tag Huffman scheme with its own
// hamming-mask-based constant-prefix detection). EncodeTags returns a
// clear error in that case rather than guessing.
package realdsrc

import (
	"fmt"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
)

// fieldCoder holds the mutable per-field encode/decode-time state a
// field's kind needs, alongside its static analyzed shape.
type fieldCoder struct {
	field   *tagField
	numeric *numericCoder // set iff field.isNumeric
	text    *textCoder    // set iff field is a free-text field (not constant, not numeric)
}

// EncodeTags analyzes and writes a block's tags to one bitstream, mirroring
// BlockCompressor::StoreTags's TagTokenizeHuffman path (the
// TagRawHuffman/mixed-formatting fallback isn't implemented).
//
// BlockCompressor::StoreTags does more than delegate to the tag encoder:
// when a block's records have varying quality lengths, it interleaves each
// record's quality length (as a fixed-width field, minQuaLength/
// maxQuaLength-derived) directly after that record's tag fields and before
// the next record's — physically part of the tag section's bitstream, not
// a separate one. qualityLengths/minQuaLength/maxQuaLength let this
// function reproduce that placement exactly; callers assembling a full
// block (see block.go) get quality lengths back from DecodeTags for
// exactly this reason.
func EncodeTags(w *bitio.Writer, tags [][]byte, qualityLengths []int, minQuaLength, maxQuaLength int) error {
	a := analyzeTags(tags)
	if a.mixedFormatting {
		return fmt.Errorf("realdsrc: tags don't share one consistent field layout (mixed formatting), which needs the unimplemented whole-tag fallback scheme")
	}

	coders, err := storeFields(w, a)
	if err != nil {
		return err
	}

	lenBits := bitLength(uint32(maxQuaLength - minQuaLength))
	isVariableLen := lenBits > 0

	for i, tag := range tags {
		fields, _ := splitTagFields(tag)
		if err := encodeNextFields(w, coders, i, fields); err != nil {
			return err
		}
		if isVariableLen {
			w.PutBits(uint32(qualityLengths[i]-minQuaLength), lenBits)
		}
	}

	w.AlignByte()
	return nil
}

// DecodeTags reconstructs numTags tags (and each one's quality length —
// see EncodeTags) from a bitstream written by EncodeTags.
func DecodeTags(r *bitio.Reader, numTags int, minQuaLength, maxQuaLength int) ([][]byte, []int, error) {
	coders, seps, err := readFields(r)
	if err != nil {
		return nil, nil, err
	}

	lenBits := bitLength(uint32(maxQuaLength - minQuaLength))
	isVariableLen := lenBits > 0

	tags := make([][]byte, numTags)
	qualityLengths := make([]int, numTags)
	for i := range tags {
		tags[i] = decodeNextFields(r, coders, seps, i)
		if isVariableLen {
			qualityLengths[i] = int(r.GetBits(lenBits)) + minQuaLength
		} else {
			qualityLengths[i] = maxQuaLength
		}
	}

	r.AlignByte()
	return tags, qualityLengths, nil
}

// storeFields mirrors TagTokenizerEncoder::StoreFields.
func storeFields(w *bitio.Writer, a *tagAnalysis) ([]*fieldCoder, error) {
	w.PutByte(byte(len(a.fields)))

	coders := make([]*fieldCoder, len(a.fields))
	for i, f := range a.fields {
		w.PutByte(f.sep)
		w.PutByte(boolBit8(f.isConstant))

		if f.isConstant {
			w.PutWord(uint32(f.len))
			for _, b := range f.data {
				w.PutByte(b)
			}
			coders[i] = &fieldCoder{field: f}
			continue
		}

		if f.isNumeric {
			nc, err := storeNumericFieldHeader(w, f)
			if err != nil {
				return nil, fmt.Errorf("field %d: %w", i, err)
			}
			coders[i] = &fieldCoder{field: f, numeric: nc}
			continue
		}

		w.PutByte(0) // is_numeric = false
		tc := storeTextFieldHeader(w, f)
		coders[i] = &fieldCoder{field: f, text: tc}
	}
	return coders, nil
}

// encodeNextFields mirrors TagTokenizerEncoder::EncodeNextFields.
func encodeNextFields(w *bitio.Writer, coders []*fieldCoder, recordIdx int, fields [][]byte) error {
	for i, fc := range coders {
		f := fc.field
		value := fields[i]

		switch {
		case f.isConstant:
			// nothing: value is implied by the header.
		case f.isNumeric:
			v, ok := isNum(value)
			if !ok {
				return fmt.Errorf("realdsrc: record %d field %d was numeric during analysis but isn't now", recordIdx, i)
			}
			encodeNumericField(w, f, fc.numeric, recordIdx, v)
		default:
			encodeTextField(w, f, fc.text, value)
		}
	}
	return nil
}

// readFields mirrors TagTokenizerDecoder::ReadFields, returning each
// field's separator alongside its coder (DecodeTags needs seps to
// reassemble a tag; EncodeTags's caller already has them from the input).
func readFields(r *bitio.Reader) ([]*fieldCoder, []byte, error) {
	n := int(r.GetByte())
	coders := make([]*fieldCoder, n)
	seps := make([]byte, n)

	for i := 0; i < n; i++ {
		sep := r.GetByte()
		seps[i] = sep
		isConstant := r.GetByte() != 0

		if isConstant {
			length := int(r.GetWord())
			data := make([]byte, length)
			for j := range data {
				data[j] = r.GetByte()
			}
			coders[i] = &fieldCoder{field: &tagField{sep: sep, isConstant: true, data: data, len: length}}
			continue
		}

		isNumeric := r.GetByte() != 0
		if isNumeric {
			f, nc, err := loadNumericFieldHeader(r)
			if err != nil {
				return nil, nil, fmt.Errorf("field %d: %w", i, err)
			}
			f.sep = sep
			coders[i] = &fieldCoder{field: f, numeric: nc}
			continue
		}

		f, tc := loadTextFieldHeader(r)
		f.sep = sep
		coders[i] = &fieldCoder{field: f, text: tc}
	}
	return coders, seps, nil
}

// decodeNextFields mirrors TagTokenizerDecoder::DecodeNextFields,
// reassembling one tag from its decoded fields and separators.
func decodeNextFields(r *bitio.Reader, coders []*fieldCoder, seps []byte, recordIdx int) []byte {
	var out []byte
	for i, fc := range coders {
		f := fc.field
		var value []byte
		switch {
		case f.isConstant:
			value = f.data
		case f.isNumeric:
			v := decodeNumericField(r, f, fc.numeric, recordIdx)
			value = []byte(itoa(v))
		default:
			value = decodeTextField(r, f, fc.text)
		}
		out = append(out, value...)
		out = append(out, seps[i])
	}
	return out[:len(out)-1] // drop the last field's separator byte (never a real one — see tag_analyze.go)
}

func itoa(v int32) string {
	if v == 0 {
		return "0"
	}
	neg := v < 0
	if neg {
		v = -v
	}
	var buf [16]byte
	pos := len(buf)
	for v > 0 {
		pos--
		buf[pos] = byte('0' + v%10)
		v /= 10
	}
	if neg {
		pos--
		buf[pos] = '-'
	}
	return string(buf[pos:])
}
