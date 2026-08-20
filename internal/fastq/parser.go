package fastq

import "errors"

// ErrNoRecords is returned when a chunk contains no valid FASTQ records.
var ErrNoRecords = errors.New("fastq: no records parsed from chunk")

// cursor walks a byte slice line by line, mirroring the Getc/Peekc/Skipc/
// SkipLine logic in dsrc::fq::FastqParser (src/FastqParser.h).
type cursor struct {
	data []byte
	pos  int
}

// skipLine advances past one line and returns its content length, excluding
// the terminator. Handles bare \n, bare \r, and \r\n exactly like upstream's
// SkipLine: it stops at the first terminator found, and returns a 0-length
// result if that terminator is the very next byte.
func (c *cursor) skipLine() int {
	n := 0
	for c.pos < len(c.data) {
		ch := c.data[c.pos]
		if ch != '\n' && ch != '\r' {
			n++
			c.pos++
			continue
		}
		if ch == '\r' && c.pos+1 < len(c.data) && c.data[c.pos+1] == '\n' {
			c.pos += 2
		} else {
			c.pos++
		}
		return n
	}
	return n
}

// readRecord ports FastqParser::ReadNextRecord.
func (c *cursor) readRecord() (Record, bool) {
	if c.pos >= len(c.data) {
		return Record{}, false
	}

	titleStart := c.pos
	titleLen := c.skipLine()
	if titleLen == 0 || c.data[titleStart] != '@' {
		return Record{}, false
	}

	seqStart := c.pos
	seqLen := c.skipLine()

	plusLen := c.skipLine()

	quaStart := c.pos
	quaLen := c.skipLine()

	if plusLen == 0 || seqLen != quaLen {
		return Record{}, false
	}

	return Record{
		Title:    c.data[titleStart : titleStart+titleLen],
		Sequence: c.data[seqStart : seqStart+seqLen],
		Quality:  c.data[quaStart : quaStart+quaLen],
	}, true
}

// ParseChunk parses every complete record out of a chunk, matching
// FastqParser::ParseFrom (src/FastqParser.cpp). It stops at the first
// malformed or incomplete record, so callers should feed it chunks that were
// produced by ChunkReader (which trims to a record boundary).
func ParseChunk(data []byte) ([]Record, StreamSizes, error) {
	c := cursor{data: data}
	var records []Record
	var sizes StreamSizes

	for c.pos < len(c.data) {
		rec, ok := c.readRecord()
		if !ok {
			break
		}
		sizes.Tag += uint64(len(rec.Title))
		sizes.Dna += uint64(len(rec.Sequence))
		sizes.Quality += uint64(len(rec.Quality))
		records = append(records, rec)
	}

	if len(records) == 0 {
		return nil, sizes, ErrNoRecords
	}
	return records, sizes, nil
}

// Analyze inspects a chunk to detect the dataset's format (color-space
// encoding, repeated '+' line, quality offset), mirroring
// FastqParser::Analyze (src/FastqParser.cpp). It returns false if the chunk
// is empty, inconsistent, or has fewer than 2 records.
func Analyze(data []byte, estimateQualityOffset bool) (DatasetType, bool) {
	c := cursor{data: data}
	var ds DatasetType
	var minQ byte = 0xFF
	var maxQ byte = 0
	recCount := 0

	for c.pos < len(c.data) {
		titleStart := c.pos
		titleLen := c.skipLine()
		if titleLen == 0 || c.data[titleStart] != '@' {
			break
		}

		seqStart := c.pos
		seqLen := c.skipLine()
		if seqLen == 0 {
			break
		}

		plusStart := c.pos
		plusRaw := c.skipLine()
		plusRep := plusRaw > 1
		if c.data[plusStart] != '+' {
			break
		}

		var quaLen int
		if estimateQualityOffset {
			quaStart := c.pos
			quaLen = c.skipLine()
			for i := 0; i < quaLen; i++ {
				q := c.data[quaStart+i]
				if q < minQ {
					minQ = q
				}
				if q > maxQ {
					maxQ = q
				}
			}
		} else {
			quaLen = c.skipLine()
			if quaLen == 0 {
				break
			}
		}
		_ = quaLen

		// Upstream reads sequence[1] unconditionally; guard the 1-byte case
		// so a short final record can't run past the slice.
		colorEnc := false
		if seqLen >= 2 {
			s1 := c.data[seqStart+1]
			colorEnc = (s1 >= '0' && s1 <= '3') || s1 == '.'
		}

		if recCount != 0 {
			if ds.ColorSpace != colorEnc {
				return ds, false
			}
			if ds.ColorSpace && c.data[seqStart] >= '0' && c.data[seqStart] <= '3' {
				return ds, false
			}
			if ds.PlusRepetition != plusRep {
				return ds, false
			}
		} else {
			ds.PlusRepetition = plusRep
			ds.ColorSpace = colorEnc
		}
		recCount++
	}

	if estimateQualityOffset {
		switch {
		case maxQ <= 74:
			if minQ >= 33 {
				ds.QualityOffset = 33 // standard Sanger / Illumina 1.8+
			}
		case maxQ <= 105:
			if minQ >= 64 {
				ds.QualityOffset = 64 // Illumina 1.3-1.8
			} else if minQ >= 59 {
				ds.QualityOffset = 59 // Solexa
			}
		}
		if ds.QualityOffset == AutoQualityOffset {
			if minQ >= 33 {
				ds.QualityOffset = 33
			} else {
				return ds, false
			}
		}
	}

	return ds, recCount > 1
}
