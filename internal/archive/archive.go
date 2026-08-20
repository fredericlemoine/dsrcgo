// Package archive is a byte-exact port of DSRC's file container format
// (src/DsrcFile.h, src/DsrcFile.cpp): a fixed 40-byte header pointing at a
// footer holding a per-block size index, followed by the concatenated
// block payloads. Writer/Reader's Create/WriteBlock/Close and
// Open/ReadNextBlock/Close mirror DsrcFileWriter's and DsrcFileReader's
// StartCompress/WriteNextChunk/FinishCompress and
// StartDecompress/ReadNextChunk/FinishDecompress.
//
// Every field, magic byte, and even the endianness quirk are verified
// against real archives produced by the original C++ dsrc binary (built
// from this repo's own src/ — see testdata/real_default_mode.dsrc, which
// is real dsrc's actual output for a small test file compressed with its
// default -m0 settings): the header uses dummy byte 0xAA and DSRC's own
// version numbers (2.0.2); the footer's block-size list is written as raw
// native-endian (little-endian on the x86/ARM machines DSRC targets)
// uint32s via a straight memcpy-style write in the C++, unlike every other
// multi-byte field in the format, which is explicitly big-endian.
//
// This container format alone is NOT enough for a real dsrc binary to
// decompress a dsrc-go archive, or vice versa: the block payloads between
// the header and footer are written by this repo's own block package,
// which is not yet byte-compatible with real DSRC's BlockCompressor
// output. This package only makes the outer envelope — magic, version,
// header/footer layout, block boundaries, and per-archive metadata —
// indistinguishable from a real one.
package archive

import (
	"encoding/binary"
	"fmt"
	"io"
	"os"
)

const (
	headerDummyByte = 0xAA
	footerDummyByte = 0xCC

	versionMajor = 2
	versionMinor = 0
	versionRev   = 2

	reservedBytes = 8
	headerSize    = 4 + reservedBytes + 3*8 + 4 // dummy+3 version bytes, footerSize, footerOffset+recordsCount+blockCount, reserved -- 40 bytes, matching DsrcFileHeader::HeaderSize

	flagPlusRepetition = 1 << 0
	flagColorSpace     = 1 << 1

	flagLossyQuality   = 1 << 0
	flagCalculateCRC32 = 1 << 1
)

// DatasetType mirrors dsrc::fq::FastqDatasetType as stored in the archive
// footer.
type DatasetType struct {
	QualityOffset  byte
	PlusRepetition bool
	ColorSpace     bool
}

// CompressionSettings mirrors the subset of dsrc::comp::CompressionSettings
// stored in the archive footer.
type CompressionSettings struct {
	DNAOrder         byte
	QualityOrder     byte
	TagPreserveFlags uint64
	Lossy            bool
	CalculateCRC32   bool
}

type header struct {
	DummyByte    uint8
	VersionMajor uint8
	VersionMinor uint8
	VersionRev   uint8
	FooterSize   uint32
	FooterOffset uint64
	RecordsCount uint64
	BlockCount   uint64
	Reserved     [reservedBytes]uint8
}

// Writer streams compressed blocks to a file, mirroring DsrcFileWriter.
type Writer struct {
	f            *os.File
	blockSizes   []uint32
	recordsCount uint64
	closed       bool

	datasetType  DatasetType
	compSettings CompressionSettings
}

// Create opens path for writing and skips past the header, whose contents
// (footer offset/size, block and record counts) aren't known until every
// block has been written.
func Create(path string) (*Writer, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	if _, err := f.Seek(headerSize, io.SeekStart); err != nil {
		f.Close()
		return nil, err
	}
	return &Writer{f: f}, nil
}

// SetDatasetType records the archive-wide FASTQ dataset metadata, mirroring
// DsrcFileWriter::SetDatasetType. Must be called before Close.
func (w *Writer) SetDatasetType(dt DatasetType) { w.datasetType = dt }

// SetCompressionSettings records the archive-wide compression settings,
// mirroring DsrcFileWriter::SetCompressionSettings. Must be called before
// Close.
func (w *Writer) SetCompressionSettings(cs CompressionSettings) { w.compSettings = cs }

// WriteBlock appends one already-compressed block to the archive.
// numRecords is only used to maintain the archive-wide record count.
func (w *Writer) WriteBlock(data []byte, numRecords int) error {
	if len(data) == 0 {
		return fmt.Errorf("archive: refusing to write an empty block")
	}
	if _, err := w.f.Write(data); err != nil {
		return err
	}
	w.blockSizes = append(w.blockSizes, uint32(len(data)))
	w.recordsCount += uint64(numRecords)
	return nil
}

// RecordsWritten returns the true cumulative record count so far. Unlike
// this value, the header field written to disk is always 0 (see Close) to
// stay byte-for-byte indistinguishable from real dsrc output — this getter
// is for callers (e.g. a CLI progress report) that want the real number.
func (w *Writer) RecordsWritten() uint64 { return w.recordsCount }

// Close writes the footer (the block size index plus dataset/compression
// metadata) and then goes back and writes the header that points at it —
// mirroring FinishCompress's write-footer-then-patch-header sequence, the
// reason Create seeks past the header instead of writing it up front.
func (w *Writer) Close() error {
	if w.closed {
		return nil
	}
	w.closed = true

	if len(w.blockSizes) == 0 {
		w.f.Close()
		return fmt.Errorf("archive: no blocks were written")
	}

	footerOffset, err := w.f.Seek(0, io.SeekCurrent)
	if err != nil {
		return err
	}
	if err := w.writeFooter(); err != nil {
		return err
	}
	footerEnd, err := w.f.Seek(0, io.SeekCurrent)
	if err != nil {
		return err
	}

	if _, err := w.f.Seek(0, io.SeekStart); err != nil {
		return err
	}
	h := header{
		DummyByte:    headerDummyByte,
		VersionMajor: versionMajor,
		VersionMinor: versionMinor,
		VersionRev:   versionRev,
		FooterSize:   uint32(footerEnd - footerOffset),
		FooterOffset: uint64(footerOffset),
		RecordsCount: 0, // matches real DSRC: DsrcFileWriter::FinishCompress never fills this in (see its own "//! TODO!!!")
		BlockCount:   uint64(len(w.blockSizes)),
	}
	for i := range h.Reserved {
		h.Reserved[i] = headerDummyByte
	}
	if err := binary.Write(w.f, binary.BigEndian, h); err != nil {
		return err
	}

	return w.f.Close()
}

func (w *Writer) writeFooter() error {
	if err := binary.Write(w.f, binary.BigEndian, uint8(footerDummyByte)); err != nil {
		return err
	}
	// Block sizes: raw native-endian uint32s (little-endian, matching
	// real DSRC's memcpy-style PutBytes((byte*)blockSizes.data(), ...) —
	// see package doc.
	if err := binary.Write(w.f, binary.LittleEndian, w.blockSizes); err != nil {
		return err
	}

	var dsFlags uint8
	if w.datasetType.ColorSpace {
		dsFlags |= flagColorSpace
	}
	if w.datasetType.PlusRepetition {
		dsFlags |= flagPlusRepetition
	}
	if err := binary.Write(w.f, binary.BigEndian, dsFlags); err != nil {
		return err
	}
	if err := binary.Write(w.f, binary.BigEndian, w.datasetType.QualityOffset); err != nil {
		return err
	}

	var csFlags uint8
	if w.compSettings.Lossy {
		csFlags |= flagLossyQuality
	}
	if w.compSettings.CalculateCRC32 {
		csFlags |= flagCalculateCRC32
	}
	if err := binary.Write(w.f, binary.BigEndian, csFlags); err != nil {
		return err
	}
	if err := binary.Write(w.f, binary.BigEndian, w.compSettings.DNAOrder); err != nil {
		return err
	}
	if err := binary.Write(w.f, binary.BigEndian, w.compSettings.QualityOrder); err != nil {
		return err
	}
	return binary.Write(w.f, binary.BigEndian, w.compSettings.TagPreserveFlags)
}

// Reader reads compressed blocks back out of a file written by Writer (or
// by real dsrc), mirroring DsrcFileReader.
type Reader struct {
	f            *os.File
	blockSizes   []uint32
	next         int
	recordsCount uint64
	datasetType  DatasetType
	compSettings CompressionSettings
}

// Open reads and validates the header and footer, then leaves the file
// positioned at the start of the first block.
func Open(path string) (*Reader, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	var h header
	if err := binary.Read(f, binary.BigEndian, &h); err != nil {
		f.Close()
		return nil, fmt.Errorf("archive: reading header: %w", err)
	}
	if h.DummyByte != headerDummyByte {
		f.Close()
		return nil, fmt.Errorf("archive: not a dsrc archive (bad magic byte)")
	}
	if h.VersionMajor != versionMajor || h.VersionMinor != versionMinor {
		f.Close()
		return nil, fmt.Errorf("archive: unsupported version %d.%d.%d", h.VersionMajor, h.VersionMinor, h.VersionRev)
	}
	if h.BlockCount == 0 {
		f.Close()
		return nil, fmt.Errorf("archive: corrupted header: zero blocks")
	}

	if _, err := f.Seek(int64(h.FooterOffset), io.SeekStart); err != nil {
		f.Close()
		return nil, err
	}
	r := &Reader{f: f, recordsCount: h.RecordsCount}
	if err := r.readFooter(int(h.BlockCount)); err != nil {
		f.Close()
		return nil, err
	}

	if _, err := f.Seek(headerSize, io.SeekStart); err != nil {
		f.Close()
		return nil, err
	}
	return r, nil
}

func (r *Reader) readFooter(blockCount int) error {
	var dummy uint8
	if err := binary.Read(r.f, binary.BigEndian, &dummy); err != nil {
		return fmt.Errorf("archive: reading footer: %w", err)
	}
	if dummy != footerDummyByte {
		return fmt.Errorf("archive: corrupted footer (bad magic byte)")
	}

	blockSizes := make([]uint32, blockCount)
	if err := binary.Read(r.f, binary.LittleEndian, &blockSizes); err != nil {
		return fmt.Errorf("archive: reading block size index: %w", err)
	}
	r.blockSizes = blockSizes

	var dsFlags uint8
	if err := binary.Read(r.f, binary.BigEndian, &dsFlags); err != nil {
		return err
	}
	r.datasetType.ColorSpace = dsFlags&flagColorSpace != 0
	r.datasetType.PlusRepetition = dsFlags&flagPlusRepetition != 0
	if err := binary.Read(r.f, binary.BigEndian, &r.datasetType.QualityOffset); err != nil {
		return err
	}

	var csFlags uint8
	if err := binary.Read(r.f, binary.BigEndian, &csFlags); err != nil {
		return err
	}
	r.compSettings.Lossy = csFlags&flagLossyQuality != 0
	r.compSettings.CalculateCRC32 = csFlags&flagCalculateCRC32 != 0
	if err := binary.Read(r.f, binary.BigEndian, &r.compSettings.DNAOrder); err != nil {
		return err
	}
	if err := binary.Read(r.f, binary.BigEndian, &r.compSettings.QualityOrder); err != nil {
		return err
	}
	return binary.Read(r.f, binary.BigEndian, &r.compSettings.TagPreserveFlags)
}

// BlockCount returns the number of blocks in the archive.
func (r *Reader) BlockCount() int { return len(r.blockSizes) }

// RecordsCount returns the header's record count. Real dsrc archives
// always report 0 here (see the package doc); archives written by this
// package's Writer do too, to stay indistinguishable.
func (r *Reader) RecordsCount() uint64 { return r.recordsCount }

// DatasetType returns the archive-wide FASTQ dataset metadata.
func (r *Reader) DatasetType() DatasetType { return r.datasetType }

// CompressionSettings returns the archive-wide compression settings.
func (r *Reader) CompressionSettings() CompressionSettings { return r.compSettings }

// ReadNextBlock returns the next block's compressed bytes, or io.EOF once
// every block has been read.
func (r *Reader) ReadNextBlock() ([]byte, error) {
	if r.next >= len(r.blockSizes) {
		return nil, io.EOF
	}
	size := r.blockSizes[r.next]
	buf := make([]byte, size)
	if _, err := io.ReadFull(r.f, buf); err != nil {
		return nil, fmt.Errorf("archive: reading block %d: %w", r.next, err)
	}
	r.next++
	return buf, nil
}

func (r *Reader) Close() error {
	return r.f.Close()
}
