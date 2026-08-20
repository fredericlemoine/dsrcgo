// Package dna ports DSRC's DNA symbol table and order-N context-modeling
// compressor (src/DnaModelerRCO.h) to Go.
package dna

// InvalidSymbol marks a byte that isn't a recognized IUPAC code.
const InvalidSymbol = 255

// AlphabetSize is the number of recognized IUPAC symbols (ACGT plus
// ambiguity/gap codes).
const AlphabetSize = 19

// iupacBases is the base->index table from RecordsProcessor::RecordsProcessor
// (src/RecordsProcessor.cpp): ACGT get indices 0-3 so the common case packs
// into 2 bits, the remaining IUPAC ambiguity codes and gap symbols follow.
var iupacBases = [AlphabetSize]byte{
	'A', 'G', 'C', 'T', 'N', 'R', 'W', 'S', 'K', 'M',
	'D', 'V', 'H', 'B', 'Y', 'X', 'U', '.', '-',
}

var toIndex [128]byte
var fromIndex [AlphabetSize]byte

func init() {
	for i := range toIndex {
		toIndex[i] = InvalidSymbol
	}
	for i, b := range iupacBases {
		toIndex[b] = byte(i)
		fromIndex[i] = b
	}
}

// ToIndex maps a FASTQ sequence byte to its symbol index. ok is false for
// any byte outside the recognized IUPAC alphabet (including lowercase).
func ToIndex(b byte) (idx byte, ok bool) {
	if b >= 128 {
		return InvalidSymbol, false
	}
	v := toIndex[b]
	return v, v != InvalidSymbol
}

// FromIndex maps a symbol index (0..AlphabetSize-1) back to its base byte.
func FromIndex(idx byte) byte {
	return fromIndex[idx]
}

// EncodeSequence maps every byte of seq to its symbol index. ok is false,
// with no partial result, if seq contains a byte outside the IUPAC alphabet.
func EncodeSequence(seq []byte) (indices []byte, ok bool) {
	out := make([]byte, len(seq))
	for i, b := range seq {
		v, valid := ToIndex(b)
		if !valid {
			return nil, false
		}
		out[i] = v
	}
	return out, true
}

// DecodeSequence maps symbol indices back to sequence bytes.
func DecodeSequence(indices []byte) []byte {
	out := make([]byte, len(indices))
	for i, v := range indices {
		out[i] = FromIndex(v)
	}
	return out
}
