// Block-level assembly, ported from BlockCompressor::Store/Read
// (src/BlockCompressor.cpp) — the piece that ties chunk metadata, tags,
// quality, and DNA together into one block, in that exact order, and runs
// the preprocessing (see preprocess.go) that couples DNA and quality.
//
// qualityOffset and plusRepetition are archive-wide settings in real DSRC
// (part of DsrcFileFooter's DatasetType, not stored per block — see
// go/internal/archive), so they're parameters here rather than something
// derived from the block itself.
package realdsrc

import (
	"bytes"
	"fmt"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/fastq"
)

// EncodeBlock compresses one chunk of raw FASTQ text (as produced by
// fastq.ChunkReader) into a real-dsrc-compatible block, mirroring
// BlockCompressor::Store.
func EncodeBlock(chunkBytes []byte, qualityOffset byte) ([]byte, error) {
	records, _, err := fastq.ParseChunk(chunkBytes)
	if err != nil {
		return nil, err
	}
	n := len(records)

	sequences := make([][]byte, n)
	qualities := make([][]byte, n)
	tags := make([][]byte, n)
	for i, r := range records {
		sequences[i] = r.Sequence
		qualities[i] = r.Quality
		tags[i] = r.Title
	}

	pre, err := PreprocessForward(sequences, qualities, qualityOffset)
	if err != nil {
		return nil, err
	}

	dnaStreams := make([][]byte, n)
	qualStreams := make([][]byte, n)
	qualityLengths := make([]int, n)
	minQ, maxQ := -1, 0
	for i, p := range pre {
		dnaStreams[i] = p.DNA
		qualStreams[i] = p.Quality
		qualityLengths[i] = len(p.Quality)
		if minQ < 0 || qualityLengths[i] < minQ {
			minQ = qualityLengths[i]
		}
		if qualityLengths[i] > maxQ {
			maxQ = qualityLengths[i]
		}
	}

	var flags uint32
	if minQ != maxQ {
		flags |= FlagVariableLength
	}

	w := bitio.NewWriter()
	WriteChunkHeader(w, ChunkHeader{
		RecordsCount: uint32(n),
		MaxQuaLength: uint32(maxQ),
		MinQuaLength: uint32(minQ),
		Flags:        flags,
		ChunkSize:    uint32(len(chunkBytes)),
	})

	if err := EncodeTags(w, tags, qualityLengths, minQ, maxQ); err != nil {
		return nil, fmt.Errorf("tags: %w", err)
	}
	if err := EncodeQuality(w, qualStreams); err != nil {
		return nil, fmt.Errorf("quality: %w", err)
	}
	if err := EncodeDNA(w, dnaStreams); err != nil {
		return nil, fmt.Errorf("dna: %w", err)
	}

	return w.Bytes(), nil
}

// DecodeBlock reconstructs raw FASTQ text from a block written by
// EncodeBlock (or by real dsrc), mirroring BlockCompressor::Read.
// plusRepetition matches the archive-wide DatasetType flag: when true,
// each record's "+" separator line repeats the tag (minus its leading
// '@') rather than being bare, mirroring BlockCompressor::ReadTags.
func DecodeBlock(data []byte, qualityOffset byte, plusRepetition bool) ([]byte, error) {
	r := bitio.NewReader(data)
	h := ReadChunkHeader(r)

	tags, qualityLengths, err := DecodeTags(r, int(h.RecordsCount), int(h.MinQuaLength), int(h.MaxQuaLength))
	if err != nil {
		return nil, fmt.Errorf("tags: %w", err)
	}

	qualStreams, err := DecodeQuality(r, qualityLengths)
	if err != nil {
		return nil, fmt.Errorf("quality: %w", err)
	}

	// A quality value >= 128 marks a smuggled ambiguous base (see
	// preprocess.go); every other position was a real DNA-stream symbol.
	dnaLengths := make([]int, len(qualStreams))
	for i, q := range qualStreams {
		for _, b := range q {
			if b < 128 {
				dnaLengths[i]++
			}
		}
	}
	dnaStreams, err := DecodeDNA(r, dnaLengths)
	if err != nil {
		return nil, fmt.Errorf("dna: %w", err)
	}

	sequences, qualities := PostprocessBackward(dnaStreams, qualStreams, qualityOffset)

	var buf bytes.Buffer
	for i := range tags {
		buf.Write(tags[i])
		buf.WriteByte('\n')
		buf.Write(sequences[i])
		buf.WriteByte('\n')
		buf.WriteByte('+')
		if plusRepetition && len(tags[i]) > 1 {
			buf.Write(tags[i][1:])
		}
		buf.WriteByte('\n')
		buf.Write(qualities[i])
		buf.WriteByte('\n')
	}
	return buf.Bytes(), nil
}
