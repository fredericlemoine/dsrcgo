// DNA encoding for real dsrc's default mode (-d0), ported from
// DnaNormalModelerProxy/DnaModelerBasicB2/DnaModelerHuffman
// (src/DnaModelerProxy.h, src/DnaModelerBasicB2.h, src/DnaModelerHuffman.cpp).
//
// DnaNormalModelerProxy::SelectSchemeId picks a scheme byte based on how
// many distinct DNA symbols (from the fixed 19-symbol IUPAC table used
// throughout this repo's dna package: A=0,G=1,C=2,T=3,N=4,...) appear
// anywhere in the block: 0 symbols -> SchemeNone, <=4 -> SchemeB2 (2-bit
// packing), otherwise -> SchemeHuffman.
//
// EncodeDNA/DecodeDNA operate on already-preprocessed raw dna-table index
// streams (see preprocess.go for why DNA can't be handled independently of
// quality) — callers run PreprocessForward/PostprocessBackward around
// these functions, not raw ASCII sequences.
//
// Encode replicates the scheme selection exactly, with one deliberate
// safety deviation: DnaModelerBasicB2::Encode calls Put2Bits(index)
// directly on the raw fixed-table index — RecordsProcessor never remaps
// DNA indices down to a dense 0..symbolCount-1 range for the B2 path the
// way it does for the Huffman path below. That means real DSRC's own
// SchemeB2 selection (symbolCount<=4) is subtly wrong whenever those <=4
// symbols aren't all below index 4 — e.g. a block containing only A, T, N,
// and one more symbol: Put2Bits would silently truncate N's index (4) to 0
// via its 2-bit mask in a release build, corrupting that base. This
// package refuses instead of reproducing that corruption: see the
// explicit index<4 check in Encode below.
package realdsrc

import (
	"fmt"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/dna"
	"github.com/fredericlemoine/dsrcgo/internal/huffman"
)

const (
	dnaSchemeB2      = 0
	dnaSchemeHuffman = 1
	dnaSchemeNone    = 255
)

// EncodeDNA writes every record's raw-index DNA stream, choosing B2 or
// Huffman to match DnaNormalModelerProxy::SelectSchemeId. streams[i] must
// contain only valid fixed-table indices (0..dna.AlphabetSize-1) — see
// PreprocessForward.
func EncodeDNA(w *bitio.Writer, streams [][]byte) error {
	var freq [dna.AlphabetSize]uint32
	for i, s := range streams {
		for _, b := range s {
			if int(b) >= dna.AlphabetSize {
				return fmt.Errorf("realdsrc: stream %d contains an invalid DNA index %d", i, b)
			}
			freq[b]++
		}
	}

	symbolCount := 0
	maxSym := -1
	for i, f := range freq {
		if f > 0 {
			symbolCount++
			maxSym = i
		}
	}

	switch {
	case symbolCount == 0:
		w.PutByte(dnaSchemeNone)
		return nil
	case symbolCount <= 4:
		if maxSym >= 4 {
			return fmt.Errorf("realdsrc: block has %d distinct DNA symbols (real dsrc's B2 threshold) but includes a non-ACGT one (fixed index %d) — refusing rather than reproducing real dsrc's silent 2-bit-truncation corruption in this case (see package doc)", symbolCount, maxSym)
		}
		return encodeDNAB2(w, streams)
	default:
		return encodeDNAHuffman(w, streams, freq)
	}
}

func encodeDNAB2(w *bitio.Writer, streams [][]byte) error {
	w.PutByte(dnaSchemeB2)
	for _, s := range streams {
		for _, b := range s {
			w.PutBits(uint32(b), 2)
		}
	}
	w.AlignByte()
	return nil
}

func encodeDNAHuffman(w *bitio.Writer, streams [][]byte, freq [dna.AlphabetSize]uint32) error {
	w.PutByte(dnaSchemeHuffman)

	// Presence bitmap over the fixed 19-symbol table, ascending raw index
	// order — mirrors the `for i in 0..MaxSymbolCount: PutBit(symbols[i] !=
	// EmptySymbol)` loop in DnaModelerHuffman::Encode. The decoder
	// reconstructs the same raw<->dense compaction from this bitmap alone.
	var rawToDense [dna.AlphabetSize]byte
	dense := 0
	for i, f := range freq {
		w.PutBit(boolBit(f > 0))
		if f > 0 {
			rawToDense[i] = byte(dense)
			dense++
		}
	}
	w.AlignByte()

	tree := huffman.NewTree()
	tree.Restart(dense)
	for _, f := range freq {
		if f > 0 {
			tree.Insert(f)
		}
	}
	codes := tree.Complete()
	storeRealHuffmanTree(w, tree)

	for _, s := range streams {
		for _, raw := range s {
			c := codes[rawToDense[raw]]
			w.PutBits(c.Code, uint(c.Len))
		}
	}
	w.AlignByte()
	return nil
}

func boolBit(b bool) uint32 {
	if b {
		return 1
	}
	return 0
}

// DecodeDNA reads len(lengths) raw-index DNA streams, one per entry of
// lengths (a stream's length may be shorter than its record's nominal read
// length — see preprocess.go).
func DecodeDNA(r *bitio.Reader, lengths []int) ([][]byte, error) {
	scheme := r.GetByte()
	switch scheme {
	case dnaSchemeNone:
		out := make([][]byte, len(lengths))
		for i, n := range lengths {
			if n != 0 {
				return nil, fmt.Errorf("realdsrc: DNA scheme is 'none' but record %d has non-zero stream length %d", i, n)
			}
		}
		return out, nil
	case dnaSchemeB2:
		return decodeDNAB2(r, lengths), nil
	case dnaSchemeHuffman:
		return decodeDNAHuffman(r, lengths)
	default:
		return nil, fmt.Errorf("realdsrc: unrecognized DNA scheme %d", scheme)
	}
}

func decodeDNAB2(r *bitio.Reader, lengths []int) [][]byte {
	out := make([][]byte, len(lengths))
	for i, n := range lengths {
		s := make([]byte, n)
		for j := range s {
			s[j] = byte(r.GetBits(2))
		}
		out[i] = s
	}
	r.AlignByte()
	return out
}

func decodeDNAHuffman(r *bitio.Reader, lengths []int) ([][]byte, error) {
	var denseToRaw []byte
	for i := 0; i < dna.AlphabetSize; i++ {
		if r.GetBit() != 0 {
			denseToRaw = append(denseToRaw, byte(i))
		}
	}
	r.AlignByte()

	tree := loadRealHuffmanTree(r)

	out := make([][]byte, len(lengths))
	for i, n := range lengths {
		s := make([]byte, n)
		for j := range s {
			dense := tree.DecodeSymbol(r)
			if int(dense) >= len(denseToRaw) {
				return nil, fmt.Errorf("realdsrc: decoded DNA symbol %d out of range (only %d symbols in this block's alphabet)", dense, len(denseToRaw))
			}
			s[j] = denseToRaw[dense]
		}
		out[i] = s
	}
	r.AlignByte()
	return out, nil
}
