// Package fastq ports DSRC's FastqParser/FastqStream layer (src/FastqParser.*,
// src/FastqStream.*) to Go. It reads a FASTQ file as fixed-size chunks aligned
// on record boundaries, then parses each chunk into records.
package fastq

// Record is a parsed FASTQ record. The slices alias the chunk they were
// parsed from — they are only valid until the chunk is reused.
//
// Mirrors dsrc::fq::FastqRecord (src/Fastq.h): the "+" separator line is
// intentionally not stored, matching upstream (only its presence/length is
// validated).
type Record struct {
	Title    []byte
	Sequence []byte
	Quality  []byte
}

// StreamSizes accumulates the per-stream byte totals produced while parsing
// a chunk, mirroring dsrc::fq::StreamsInfo (src/Common.h).
type StreamSizes struct {
	Tag     uint64
	Dna     uint64
	Quality uint64
}

// DatasetType mirrors dsrc::fq::FastqDatasetType (src/Common.h).
type DatasetType struct {
	QualityOffset  uint32
	PlusRepetition bool
	ColorSpace     bool
}

const (
	AutoQualityOffset    uint32 = 0
	DefaultQualityOffset uint32 = 33
)
