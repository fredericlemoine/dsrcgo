// Quality's Truncated scheme, ported from QualityPositionModelerTruncated
// (src/QualityPositionModeler.cpp) — used when a block's quality strings
// tend to end in a trailing run of Phred-2 ("poor/undetermined") calls: the
// trailing run is stripped before building per-position statistics and
// codes, and only the truncated prefix is actually Huffman-coded per
// record; the stripped tail is reconstructed on decode by refilling with
// the sentinel value 2 rather than storing it at all.
package realdsrc

import (
	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/huffman"
)

// bitLength mirrors core::bit_length (src/utils.h): the number of bits
// needed to represent x, i.e. the smallest n with x < 2^n (so
// bitLength(0) == 0, unlike this package's bitsForCount, which returns the
// bits needed to distinguish n *values* 0..n-1).
func bitLength(x uint32) uint {
	for i := uint(0); i < 32; i++ {
		if x < (uint32(1) << i) {
			return i
		}
	}
	return 64
}

// truncatedLen mirrors the truncatedLen computation in
// LosslessRecordsProcessor::ProcessForward: one past the index of the last
// quality value that isn't the sentinel 2, i.e. q with any trailing run of
// 2s stripped (interior 2s don't matter, only a trailing run does).
func truncatedLen(q []byte) int {
	if len(q) == 0 {
		return 0
	}
	last := 0
	for i, b := range q {
		if b != hashSymbolNormalQ {
			last = i
		}
	}
	return last + 1
}

func encodeQualityTruncated(w *bitio.Writer, qualityStreams [][]byte, st qualityStats) error {
	maxLen := st.maxLength

	rawToDense, dense := qualityRawToDense(st.freq)

	tLens := make([]int, len(qualityStreams))
	for i, q := range qualityStreams {
		tLens[i] = truncatedLen(q)
	}

	posFreq := make([][]uint32, maxLen)
	for j := range posFreq {
		posFreq[j] = make([]uint32, dense)
	}
	for i, q := range qualityStreams {
		for j := 0; j < tLens[i]; j++ {
			posFreq[j][rawToDense[q[j]]]++
		}
	}

	w.PutWord(uint32(maxLen))
	storeQualitySymbolBitmap(w, st.freq)

	trees, allCodes := buildPositionTrees(posFreq, dense)
	for _, tr := range trees {
		storeRealHuffmanTree(w, tr)
	}

	variableLength := st.minLength != st.maxLength
	maxBitLength := bitLength(uint32(maxLen))
	w.PutBit(boolBit(variableLength))

	for i, q := range qualityStreams {
		qLen, tLen := len(q), tLens[i]
		w.PutBit(boolBit(qLen != tLen))
		if qLen != tLen {
			bitLen := maxBitLength
			if variableLength {
				bitLen = bitLength(uint32(qLen))
			}
			w.PutBits(uint32(tLen), bitLen)
		}
		for j := 0; j < tLen; j++ {
			c := allCodes[j][rawToDense[q[j]]]
			w.PutBits(c.Code, uint(c.Len))
		}
	}
	w.AlignByte()
	return nil
}

func decodeQualityTruncated(r *bitio.Reader, lengths []int) ([][]byte, error) {
	maxLen := int(r.GetWord())
	denseToRaw := loadQualitySymbolBitmap(r)

	trees := make([]*realHuffmanTree, maxLen)
	for j := 0; j < maxLen; j++ {
		trees[j] = loadRealHuffmanTree(r)
	}

	variableLength := r.GetBit() != 0
	maxBitLength := bitLength(uint32(maxLen))

	out := make([][]byte, len(lengths))
	for i, qLen := range lengths {
		thLen := qLen
		if r.GetBit() != 0 {
			bitLen := maxBitLength
			if variableLength {
				bitLen = bitLength(uint32(qLen))
			}
			thLen = int(r.GetBits(bitLen))
		}

		q := make([]byte, qLen)
		for j := 0; j < thLen; j++ {
			dense := trees[j].DecodeSymbol(r)
			if int(dense) >= len(denseToRaw) {
				return nil, errQualitySymbolRange(dense, j, len(denseToRaw))
			}
			q[j] = denseToRaw[dense]
		}
		for j := thLen; j < qLen; j++ {
			q[j] = hashSymbolNormalQ
		}
		out[i] = q
	}
	r.AlignByte()
	return out, nil
}

// qualityRawToDense builds the raw-quality-value -> dense-index table in
// ascending value order, matching QualityStats.symbols's compaction.
func qualityRawToDense(freq [maxQualitySymbolCount]uint32) (table [maxQualitySymbolCount]byte, dense int) {
	for v, f := range freq {
		if f > 0 {
			table[v] = byte(dense)
			dense++
		}
	}
	return table, dense
}

func storeQualitySymbolBitmap(w *bitio.Writer, freq [maxQualitySymbolCount]uint32) {
	for v := 0; v < maxQualitySymbolCount; v++ {
		w.PutBit(boolBit(freq[v] > 0))
	}
}

func loadQualitySymbolBitmap(r *bitio.Reader) []byte {
	var denseToRaw []byte
	for v := 0; v < maxQualitySymbolCount; v++ {
		if r.GetBit() != 0 {
			denseToRaw = append(denseToRaw, byte(v))
		}
	}
	return denseToRaw
}

// buildPositionTrees builds one real-format Huffman tree per position from
// per-position symbol frequencies, shared by the Plain and Truncated
// schemes.
func buildPositionTrees(posFreq [][]uint32, dense int) ([]*huffman.Tree, [][]huffman.Code) {
	trees := make([]*huffman.Tree, len(posFreq))
	codes := make([][]huffman.Code, len(posFreq))
	for j, freqs := range posFreq {
		tr := huffman.NewTree()
		tr.Restart(dense)
		for _, f := range freqs {
			tr.Insert(f)
		}
		codes[j] = tr.Complete()
		trees[j] = tr
	}
	return trees, codes
}
