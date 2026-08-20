// Quality encoding for real dsrc's default mode (-q0), ported from
// QualityNormalModelerProxy and QualityPositionModelerPlain
// (src/QualityModelerProxy.h, src/QualityPositionModeler.cpp).
//
// QualityNormalModelerProxy::SelectSchemeId picks among three schemes
// using aggregate stats gathered while preprocessing every record in a
// block: QualityRle when thLength/rleLength > 1.25, QualityTruncated when
// rawLength/thLength > 1.10, else QualityPlain. EncodeQuality replicates
// that selection exactly (via computeQualityStats/selectQualityScheme) but
// only implements the Plain scheme's actual encoding — Truncated and RLE
// return a clear error instead of silently producing incompatible output.
// Plain is upstream's fallback case, so it's the common one for data
// without either long low-complexity quality tails (which favor
// Truncated) or long runs of a repeated quality value (which favor RLE).
//
// Plain itself: one independent Huffman tree per read offset (real
// DSRC's actual dual-representation tree format — see realhuffman.go —
// not go/internal/quality's PositionModel, which deliberately uses a
// different, simpler raw<->dense table; see that package's doc comment
// for why), built from the distribution of quality values observed at
// that offset across the whole block.
package realdsrc

import (
	"fmt"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
)

const (
	qualitySchemePlain     = 0
	qualitySchemeTruncated = 1
	qualitySchemeRLE       = 2

	maxQualitySymbolCount = 256
	hashSymbolNormalQ     = 2 // IRecordsProcessor::HashSymbolNormal, used by the truncation heuristic below
)

// qualityStats mirrors the subset of QualityStats needed for scheme
// selection: rawLength, thLength, rleLength (see the package doc for the
// selection formula), plus per-block symbol frequencies and min/max
// length needed by whichever scheme actually gets used.
type qualityStats struct {
	rawLength, thLength, rleLength uint64
	minLength, maxLength           int
	freq                           [maxQualitySymbolCount]uint32
}

// computeQualityStats mirrors the per-record loop in
// LosslessRecordsProcessor::ProcessForward that accumulates
// QualityStats — everything except the DNA/smuggling side of that
// function, which preprocess.go already handled to produce these streams.
func computeQualityStats(qualityStreams [][]byte) qualityStats {
	var st qualityStats
	st.minLength = -1

	for _, q := range qualityStreams {
		var prevSym int = -1 // EmptySymbol sentinel (255) can't collide with byte values here since we just need "no previous symbol yet"
		curQThLen := 0

		for i, b := range q {
			st.freq[b]++

			if int(b) != prevSym {
				st.rleLength++
			}
			if b != hashSymbolNormalQ {
				curQThLen = i
			}
			prevSym = int(b)
		}

		if prevSym == hashSymbolNormalQ && st.rleLength > 0 {
			st.rleLength--
		}

		st.rawLength += uint64(len(q))
		st.thLength += uint64(curQThLen)

		if st.minLength < 0 || len(q) < st.minLength {
			st.minLength = len(q)
		}
		if len(q) > st.maxLength {
			st.maxLength = len(q)
		}
	}
	return st
}

// selectQualityScheme mirrors QualityNormalModelerProxy::SelectSchemeId.
func selectQualityScheme(st qualityStats) byte {
	if st.rleLength > 0 && float64(st.thLength)/float64(st.rleLength) > 1.25 {
		return qualitySchemeRLE
	}
	if st.thLength > 0 && float64(st.rawLength)/float64(st.thLength) > 1.10 {
		return qualitySchemeTruncated
	}
	return qualitySchemePlain
}

// EncodeQuality writes every record's quality stream (see preprocess.go —
// this is the post-smuggling stream, same length as the nominal quality
// string), choosing a scheme the same way real dsrc does.
func EncodeQuality(w *bitio.Writer, qualityStreams [][]byte) error {
	st := computeQualityStats(qualityStreams)
	scheme := selectQualityScheme(st)

	w.PutByte(scheme)
	switch scheme {
	case qualitySchemePlain:
		return encodeQualityPlain(w, qualityStreams, st)
	case qualitySchemeTruncated:
		return encodeQualityTruncated(w, qualityStreams, st)
	default:
		return encodeQualityRLE(w, qualityStreams)
	}
}

func encodeQualityPlain(w *bitio.Writer, qualityStreams [][]byte, st qualityStats) error {
	maxLen := st.maxLength
	rawToDense, dense := qualityRawToDense(st.freq)

	// Per-position symbol frequency, used to build one tree per offset.
	posFreq := make([][]uint32, maxLen)
	for j := range posFreq {
		posFreq[j] = make([]uint32, dense)
	}
	for _, q := range qualityStreams {
		for j, b := range q {
			posFreq[j][rawToDense[b]]++
		}
	}

	w.PutWord(uint32(maxLen)) // StoreStatsData
	storeQualitySymbolBitmap(w, st.freq)

	trees, allCodes := buildPositionTrees(posFreq, dense)
	for _, tr := range trees {
		storeRealHuffmanTree(w, tr)
	}

	// EncodeRecords: per record, per position, the position's Huffman code
	// for that quality value's dense index.
	for _, q := range qualityStreams {
		for j, b := range q {
			c := allCodes[j][rawToDense[b]]
			w.PutBits(c.Code, uint(c.Len))
		}
	}
	w.AlignByte()
	return nil
}

// DecodeQuality reads a quality section written by EncodeQuality (or by
// real dsrc, for the Plain scheme) and returns one quality stream per
// entry of lengths.
func DecodeQuality(r *bitio.Reader, lengths []int) ([][]byte, error) {
	scheme := r.GetByte()
	switch scheme {
	case qualitySchemePlain:
		return decodeQualityPlain(r, lengths)
	case qualitySchemeTruncated:
		return decodeQualityTruncated(r, lengths)
	case qualitySchemeRLE:
		return decodeQualityRLE(r, lengths)
	default:
		return nil, fmt.Errorf("realdsrc: unrecognized quality scheme %d", scheme)
	}
}

func decodeQualityPlain(r *bitio.Reader, lengths []int) ([][]byte, error) {
	maxLen := int(r.GetWord())
	denseToRaw := loadQualitySymbolBitmap(r)

	trees := make([]*realHuffmanTree, maxLen)
	for j := 0; j < maxLen; j++ {
		trees[j] = loadRealHuffmanTree(r)
	}

	out := make([][]byte, len(lengths))
	for i, n := range lengths {
		q := make([]byte, n)
		for j := 0; j < n; j++ {
			dense := trees[j].DecodeSymbol(r)
			if int(dense) >= len(denseToRaw) {
				return nil, errQualitySymbolRange(dense, j, len(denseToRaw))
			}
			q[j] = denseToRaw[dense]
		}
		out[i] = q
	}
	r.AlignByte()
	return out, nil
}

func errQualitySymbolRange(dense int32, position, alphabetSize int) error {
	return fmt.Errorf("realdsrc: decoded quality symbol %d out of range at position %d (alphabet has %d symbols)", dense, position, alphabetSize)
}
