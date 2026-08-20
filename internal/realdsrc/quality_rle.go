// Quality's RLE scheme, ported from QualityRLEModeler
// (src/QualityRLEModeler.cpp) — used when a block's quality values form
// long runs of a repeated symbol. Every record's quality stream is
// concatenated into one continuous run-length sequence (record boundaries
// don't reset a run — the stream really is treated block-wide), each run
// capped at length 255 (encoded as length-1, 0..254). Each run's symbol is
// then coded against a Huffman tree selected by the *previous* run's
// symbol, and each run's length against a tree selected by its *own*
// symbol — a first-order context model over runs, not raw positions.
//
// A degenerate case gets its own non-Huffman path: if the whole block has
// only one distinct quality value, there's nothing to condition on, so
// EncodeRuns instead relies on an invariant upstream asserts rather than
// proves in general — a single repeated symbol split only by the 255-run
// cap produces at most 2 distinct run lengths (the capped length, repeated,
// and one shorter remainder) — and stores just the first run's length,
// inferring the rest structurally. This package checks that invariant
// explicitly rather than assuming it (see decodeQualityRLESingleSymbol).
package realdsrc

import (
	"fmt"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/huffman"
)

// maxRunLength mirrors QualityRLEModeler::MaxLengthSymbol: a run's encoded
// length (actual run length - 1) is capped at this value, so the longest
// single run is 255 quality values.
const maxRunLength = 254

type qualityRun struct {
	symbol byte // raw quality value
	length byte // actual run length - 1
}

// computeRuns run-length-encodes every record's quality stream as one
// continuous sequence, mirroring QualityRLEModeler::EncodeRecords. Callers
// must pass at least one non-empty stream (always true for valid FASTQ
// data, matching every other assumption in this package).
func computeRuns(qualityStreams [][]byte) []qualityRun {
	var runs []qualityRun
	prevSym := -1
	curLen := 0

	for _, q := range qualityStreams {
		for _, b := range q {
			if int(b) == prevSym && curLen < maxRunLength {
				curLen++
				continue
			}
			if prevSym != -1 {
				runs = append(runs, qualityRun{byte(prevSym), byte(curLen)})
			}
			curLen = 0
			prevSym = int(b)
		}
	}
	runs = append(runs, qualityRun{byte(prevSym), byte(curLen)})
	return runs
}

func encodeQualityRLE(w *bitio.Writer, qualityStreams [][]byte) error {
	runs := computeRuns(qualityStreams)

	var qF, lF [maxQualitySymbolCount]uint32
	for _, r := range runs {
		qF[r.symbol]++
		lF[r.length]++
	}
	qRawToDense, qDense := qualityRawToDense(qF)
	lRawToDense, lDense := qualityRawToDense(lF)

	w.PutWord(uint32(len(runs))) // StoreStatsData
	storeQualitySymbolBitmap(w, qF)
	storeQualitySymbolBitmap(w, lF)

	if qDense == 1 {
		if lDense > 2 {
			return fmt.Errorf("realdsrc: RLE single-symbol block has %d distinct run lengths, expected <=2 (real dsrc's own degenerate-case encoding assumes this)", lDense)
		}
		if lDense > 1 {
			w.AlignByte()
			w.PutByte(lRawToDense[runs[0].length])
		}
		w.AlignByte()
		return nil
	}

	qContexts, lContexts, qCodes, lCodes := buildRunContexts(runs, qRawToDense, qDense, lRawToDense, lDense)
	for i := 0; i < qDense; i++ {
		storeRealHuffmanTree(w, qContexts[i])
		storeRealHuffmanTree(w, lContexts[i])
	}

	prev := 0
	for _, r := range runs {
		q := qRawToDense[r.symbol]
		l := lRawToDense[r.length]
		qc := qCodes[prev][q]
		lc := lCodes[q][l]
		w.PutBits(qc.Code, uint(qc.Len))
		w.PutBits(lc.Code, uint(lc.Len))
		prev = int(q)
	}
	w.AlignByte()
	return nil
}

// buildRunContexts builds one Huffman tree per distinct previous-symbol
// context for run symbols, and one per distinct symbol context for run
// lengths, mirroring QualityRLEModeler::ComputeHuffmanContext. The very
// first run is coded under context 0 regardless of its actual symbol
// (`uchar prev = 0;` in upstream) — a fixed, not a sentinel, starting
// context, so encoder and decoder simply agree on it structurally.
func buildRunContexts(runs []qualityRun, qRawToDense [maxQualitySymbolCount]byte, qDense int, lRawToDense [maxQualitySymbolCount]byte, lDense int) (qContexts, lContexts []*huffman.Tree, qCodes, lCodes [][]huffman.Code) {
	qF := make([][]uint32, qDense)
	lF := make([][]uint32, qDense)
	for i := range qF {
		qF[i] = make([]uint32, qDense)
		lF[i] = make([]uint32, lDense)
	}

	prev := 0
	for _, r := range runs {
		q := qRawToDense[r.symbol]
		l := lRawToDense[r.length]
		qF[prev][q]++
		lF[q][l]++
		prev = int(q)
	}

	qContexts = make([]*huffman.Tree, qDense)
	lContexts = make([]*huffman.Tree, qDense)
	qCodes = make([][]huffman.Code, qDense)
	lCodes = make([][]huffman.Code, qDense)
	for i := 0; i < qDense; i++ {
		qt := huffman.NewTree()
		qt.Restart(qDense)
		for _, f := range qF[i] {
			qt.Insert(f)
		}
		qCodes[i] = qt.Complete()
		qContexts[i] = qt

		lt := huffman.NewTree()
		lt.Restart(lDense)
		for _, f := range lF[i] {
			lt.Insert(f)
		}
		lCodes[i] = lt.Complete()
		lContexts[i] = lt
	}
	return qContexts, lContexts, qCodes, lCodes
}

func decodeQualityRLE(r *bitio.Reader, lengths []int) ([][]byte, error) {
	runLength := int(r.GetWord())
	qDenseToRaw := loadQualitySymbolBitmap(r)
	lDenseToRaw := loadQualitySymbolBitmap(r)
	qDense, lDense := len(qDenseToRaw), len(lDenseToRaw)

	if qDense == 0 || lDense == 0 {
		return nil, fmt.Errorf("realdsrc: RLE quality section has an empty symbol alphabet")
	}

	var runs []qualityRun
	var err error
	if qDense == 1 {
		runs, err = decodeQualityRLESingleSymbol(r, runLength, qDenseToRaw, lDenseToRaw)
	} else {
		runs, err = decodeQualityRLEContexts(r, runLength, qDenseToRaw, lDenseToRaw, qDense)
	}
	if err != nil {
		return nil, err
	}
	r.AlignByte()

	flat := make([]byte, 0, runLength*2)
	for _, run := range runs {
		n := int(run.length) + 1
		for k := 0; k < n; k++ {
			flat = append(flat, run.symbol)
		}
	}

	out := make([][]byte, len(lengths))
	off := 0
	for i, n := range lengths {
		if off+n > len(flat) {
			return nil, fmt.Errorf("realdsrc: RLE run data (%d symbols) is shorter than the requested %d quality bytes", len(flat), off+n)
		}
		out[i] = flat[off : off+n]
		off += n
	}
	return out, nil
}

func decodeQualityRLESingleSymbol(r *bitio.Reader, runLength int, qDenseToRaw, lDenseToRaw []byte) ([]qualityRun, error) {
	if len(lDenseToRaw) > 2 {
		return nil, fmt.Errorf("realdsrc: RLE single-symbol block has %d distinct run lengths, expected <=2", len(lDenseToRaw))
	}

	qSym := qDenseToRaw[0]
	var lBegin, lEnd byte
	if len(lDenseToRaw) > 1 {
		r.AlignByte()
		lBegin = lDenseToRaw[r.GetByte()]
		lEnd = lDenseToRaw[0]
		if lEnd == lBegin {
			lEnd = lDenseToRaw[1]
		}
	} else {
		lBegin = lDenseToRaw[0]
		lEnd = lBegin
	}

	runs := make([]qualityRun, runLength)
	for i := range runs {
		runs[i] = qualityRun{qSym, lBegin}
	}
	if runLength > 0 {
		runs[runLength-1].length = lEnd
	}
	return runs, nil
}

func decodeQualityRLEContexts(r *bitio.Reader, runLength int, qDenseToRaw, lDenseToRaw []byte, qDense int) ([]qualityRun, error) {
	qContexts := make([]*realHuffmanTree, qDense)
	lContexts := make([]*realHuffmanTree, qDense)
	for i := 0; i < qDense; i++ {
		qContexts[i] = loadRealHuffmanTree(r)
		lContexts[i] = loadRealHuffmanTree(r)
	}

	runs := make([]qualityRun, runLength)
	prev := 0
	for i := range runs {
		qd := qContexts[prev].DecodeSymbol(r)
		if int(qd) >= len(qDenseToRaw) {
			return nil, errQualitySymbolRange(qd, i, len(qDenseToRaw))
		}
		ld := lContexts[qd].DecodeSymbol(r)
		if int(ld) >= len(lDenseToRaw) {
			return nil, errQualitySymbolRange(ld, i, len(lDenseToRaw))
		}
		runs[i] = qualityRun{qDenseToRaw[qd], lDenseToRaw[ld]}
		prev = int(qd)
	}
	return runs, nil
}
