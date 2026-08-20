package fastq

import "io"

// DefaultChunkSize mirrors dsrc::core::DataChunk::DefaultBufferSize (src/Buffer.h).
const DefaultChunkSize = 1 << 20

// tailScanWindow mirrors IFastqStreamReader::SwapBufferSize (src/FastqStream.h):
// the trailing region of a full buffer that gets re-scanned for a clean
// record boundary and carried over into the next chunk.
const tailScanWindow = 1 << 13

// ChunkReader splits a FASTQ stream into chunks that always end on a record
// boundary, so each chunk can be parsed independently. It ports
// IFastqStreamReader::ReadNextChunk (src/FastqStream.cpp).
type ChunkReader struct {
	r         io.Reader
	chunkSize int
	carry     []byte
	carryLen  int
	eof       bool
	usesCRLF  bool
}

// NewChunkReader wraps r. chunkSize <= 0 uses DefaultChunkSize; it must be
// larger than tailScanWindow.
func NewChunkReader(r io.Reader, chunkSize int) *ChunkReader {
	if chunkSize <= 0 {
		chunkSize = DefaultChunkSize
	}
	return &ChunkReader{
		r:         r,
		chunkSize: chunkSize,
		carry:     make([]byte, tailScanWindow),
	}
}

// ReadNextChunk returns the next record-aligned chunk, or io.EOF once the
// stream is exhausted. The returned slice is only valid until the next call.
func (cr *ChunkReader) ReadNextChunk() ([]byte, error) {
	if cr.eof {
		return nil, io.EOF
	}

	buf := make([]byte, cr.chunkSize)
	size := 0
	if cr.carryLen > 0 {
		size = copy(buf, cr.carry[:cr.carryLen])
		cr.carryLen = 0
	}

	n, err := io.ReadFull(cr.r, buf[size:])
	total := size + n

	switch err {
	case nil:
		// Buffer filled completely; more data may follow. Trim back to the
		// last clean record boundary and carry the remainder forward.
		boundary := cr.nextRecordBoundary(buf, len(buf)-tailScanWindow)
		chunkLen := boundary - 1
		if cr.usesCRLF {
			chunkLen--
		}
		cr.carryLen = copy(cr.carry, buf[boundary:])
		return buf[:chunkLen], nil

	case io.EOF, io.ErrUnexpectedEOF:
		cr.eof = true
		chunkLen := total - 1 // drop the final line terminator byte
		if cr.usesCRLF {
			chunkLen--
		}
		if chunkLen <= 0 {
			return nil, io.EOF
		}
		return buf[:chunkLen], nil

	default:
		return nil, err
	}
}

// skipToEOL advances pos to the line terminator (or end of buffer), mirroring
// IFastqStreamReader::SkipToEol.
func skipToEOL(data []byte, pos int) (next int, crlf bool) {
	size := len(data)
	for pos < size && data[pos] != '\n' && data[pos] != '\r' {
		pos++
	}
	if pos < size && data[pos] == '\r' && pos+1 < size && data[pos+1] == '\n' {
		pos++
		crlf = true
	}
	return pos, crlf
}

// nextRecordBoundary finds the start of the next record at or after pos,
// disambiguating a '@' that begins a quality line from one that begins a
// title line by checking whether the line two rows down starts with '+'.
// Mirrors IFastqStreamReader::GetNextRecordPos.
func (cr *ChunkReader) nextRecordBoundary(data []byte, pos int) int {
	size := len(data)
	adv := func(p int) int {
		next, crlf := skipToEOL(data, p)
		if crlf {
			cr.usesCRLF = true
		}
		return next + 1
	}

	pos = adv(pos)
	for pos < size && data[pos] != '@' {
		pos = adv(pos)
	}
	if pos >= size {
		return size
	}
	pos0 := pos

	pos = adv(pos)
	if pos >= size {
		return size
	}
	if data[pos] == '@' {
		return pos // pos0 was a quality line that happened to start with '@'
	}

	pos = adv(pos)
	if pos >= size {
		return size
	}
	// data[pos] should be '+', confirming pos0 was a title line.
	return pos0
}
