// Package realdsrc is a byte-exact port of DSRC's default-mode
// (-m0 = -d0 -q0, no -l, no -c, no -f) block format
// (src/BlockCompressor.cpp), verified continuously against archives
// produced by the real C++ dsrc binary — built from this repo's own src/
// (see testdata/) — rather than derived from reading the source alone.
//
// This is a from-scratch implementation, not a reuse of the block/dna/
// quality/tag packages built earlier in this port: those were built to be
// algorithmically faithful but were explicitly NOT byte-compatible with
// real DSRC (documented in each of their package comments). Achieving
// real compatibility means matching upstream's exact bit layout, which in
// several places differs from the cleaner approach those packages took —
// e.g. real DSRC's Huffman coder uses two different in-memory tree shapes
// between encode and decode sides, where go/internal/huffman deliberately
// uses one.
//
// Status: chunk metadata and DNA's B2 (2-bit ACGT packing) scheme are
// implemented and verified bit-exact against real dsrc output. Quality
// and tag encoding, and DNA's Huffman scheme (for sequences with IUPAC
// ambiguity codes), are not yet implemented — see the doc comments on
// each file for what's covered.
package realdsrc

import "github.com/fredericlemoine/dsrcgo/internal/bitio"

// bitsForCount returns the smallest b with 2^b >= n, matching upstream's
// int_log-based bits_per_id formula (core::int_log + a power-of-2 bump) —
// see realhuffman.go for the derivation showing these are equivalent.
func bitsForCount(n int) uint {
	b := uint(0)
	for (1 << b) < n {
		b++
	}
	return b
}

// Chunk header flag bits, mirroring BlockCompressor::FastqBlockFlags
// (src/BlockCompressor.h).
const (
	FlagDeltaConstant        = 1 << 0
	FlagVariableLength       = 1 << 1
	FlagMixedFieldFormatting = 1 << 2
)

// ChunkHeader mirrors the subset of dsrc::comp::ChunkHeader this package
// implements: no color-space (FlagDeltaConstant) and no CRC32 support yet.
type ChunkHeader struct {
	RecordsCount uint32
	MaxQuaLength uint32
	MinQuaLength uint32
	Flags        uint32
	ChunkSize    uint32
}

// WriteChunkHeader mirrors BlockCompressor::StoreMetaData for a dataset
// that is not color-space and an archive with CRC32 disabled.
func WriteChunkHeader(w *bitio.Writer, h ChunkHeader) {
	w.PutWord(h.RecordsCount)
	w.PutWord(h.MaxQuaLength)
	w.PutWord(h.Flags)
	w.PutWord(h.ChunkSize)
	if h.Flags&FlagVariableLength != 0 {
		w.PutWord(h.MinQuaLength)
	}
	w.AlignByte() // FlushPartialWordBuffer; a no-op here since every write above is already word-aligned
}

// ReadChunkHeader mirrors BlockCompressor::ReadMetaData for the same
// restricted case as WriteChunkHeader.
func ReadChunkHeader(r *bitio.Reader) ChunkHeader {
	var h ChunkHeader
	h.RecordsCount = r.GetWord()
	h.MaxQuaLength = r.GetWord()
	h.Flags = r.GetWord()
	h.ChunkSize = r.GetWord()
	if h.Flags&FlagVariableLength != 0 {
		h.MinQuaLength = r.GetWord()
	} else {
		h.MinQuaLength = h.MaxQuaLength
	}
	r.AlignByte()
	return h
}
