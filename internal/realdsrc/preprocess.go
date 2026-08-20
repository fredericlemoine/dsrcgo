// Package-level preprocessing shared by DNA and quality encoding, ported
// from LosslessRecordsProcessor::ProcessForward/ProcessBackward
// (src/RecordsProcessor.cpp) — the piece BlockCompressor::PreprocessRecords/
// PostprocessRecords runs before/after StoreDNA/StoreQuality.
//
// DNA and quality are NOT independent streams in real DSRC's lossless
// mode: whenever a non-ACGT base's offset-adjusted quality score is < 7,
// that base is dropped from the DNA stream entirely and "smuggled" into
// the quality stream instead, as a value >= 128 encoding both the base's
// identity and its (3-bit) quality score. This is why a record's DNA
// stream can be shorter than its nominal read length, and why decoding
// must read quality before DNA and then recombine them — matching
// BlockCompressor::Read's actual section order (tags, quality, DNA).
package realdsrc

import (
	"fmt"

	"github.com/fredericlemoine/dsrcgo/internal/dna"
)

// hashSymbolNormal mirrors IRecordsProcessor::HashSymbolNormal — a
// (offset-adjusted) quality value reserved as a sentinel elsewhere in
// upstream's truncation logic; not otherwise used by this file.
const hashSymbolNormal = 2

// PreprocessedRecord holds one record after the forward pass: DNA indices
// (0..18, fixed dna-table values) for bases that stayed in the DNA stream,
// and quality values (one per original position, 0-based post-offset,
// possibly >= 128 for a smuggled ambiguous base) for the quality stream.
type PreprocessedRecord struct {
	DNA     []byte
	Quality []byte
}

// PreprocessForward runs the forward pass over every record in a block.
// qualityOffset is the dataset's Phred offset (33 for Sanger/Illumina
// 1.8+, the value real dsrc auto-detects for essentially all modern data).
func PreprocessForward(sequences, qualities [][]byte, qualityOffset byte) ([]PreprocessedRecord, error) {
	out := make([]PreprocessedRecord, len(sequences))
	for i := range sequences {
		seq, qual := sequences[i], qualities[i]
		if len(seq) != len(qual) {
			return nil, fmt.Errorf("realdsrc: record %d has sequence length %d but quality length %d", i, len(seq), len(qual))
		}

		dnaOut := make([]byte, 0, len(seq))
		qualOut := make([]byte, len(qual))

		for j := range seq {
			rawIdx, ok := dna.ToIndex(seq[j])
			if !ok {
				return nil, fmt.Errorf("realdsrc: record %d contains a byte outside the IUPAC alphabet at position %d", i, j)
			}
			if qual[j] < qualityOffset {
				return nil, fmt.Errorf("realdsrc: record %d quality byte 0x%02x at position %d is below the quality offset %d", i, qual[j], j, qualityOffset)
			}
			qRaw := int(qual[j]) - int(qualityOffset)

			if rawIdx > 3 && qRaw < 7 {
				qualOut[j] = byte(128 + ((int(rawIdx) - 3 + 1) << 3) - 16 + qRaw)
			} else {
				dnaOut = append(dnaOut, rawIdx)
				qualOut[j] = byte(qRaw)
			}
		}

		out[i] = PreprocessedRecord{DNA: dnaOut, Quality: qualOut}
	}
	return out, nil
}

// PostprocessBackward recombines decoded DNA and quality streams back into
// sequence/quality byte strings, mirroring
// LosslessRecordsProcessor::ProcessBackward. dnaStream and qualityStream
// must correspond 1:1 (same record order) with qualityStream's lengths
// matching each record's true (nominal) length.
func PostprocessBackward(dnaStream, qualityStream [][]byte, qualityOffset byte) (sequences, qualities [][]byte) {
	sequences = make([][]byte, len(qualityStream))
	qualities = make([][]byte, len(qualityStream))

	for i, qual := range qualityStream {
		seq := dnaStream[i]
		seqi := len(seq) - 1

		outSeq := make([]byte, len(qual))
		outQual := make([]byte, len(qual))

		for j := len(qual) - 1; j >= 0; j-- {
			qval := int(qual[j])
			var seqval byte
			if qval >= 128 {
				seqval = byte((qval-128+16)/8 + 3 - 1)
				qval &= 7
			} else {
				seqval = seq[seqi]
				seqi--
			}
			outSeq[j] = dna.FromIndex(seqval)
			outQual[j] = qualityOffset + byte(qval)
		}

		sequences[i] = outSeq
		qualities[i] = outQual
	}
	return sequences, qualities
}
