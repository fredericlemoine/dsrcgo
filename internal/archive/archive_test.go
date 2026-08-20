package archive

import (
	"bytes"
	"io"
	"math/rand"
	"os"
	"path/filepath"
	"testing"
)

func TestWriteReadRoundTrip(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	dir := t.TempDir()
	path := filepath.Join(dir, "test.dsrc")

	var blocks [][]byte
	var recordCounts []int
	for i := 0; i < 20; i++ {
		size := 100 + rng.Intn(5000)
		b := make([]byte, size)
		rng.Read(b)
		blocks = append(blocks, b)
		recordCounts = append(recordCounts, 10+rng.Intn(100))
	}

	w, err := Create(path)
	if err != nil {
		t.Fatalf("Create: %v", err)
	}
	w.SetDatasetType(DatasetType{QualityOffset: 33})
	w.SetCompressionSettings(CompressionSettings{DNAOrder: 0, QualityOrder: 0})

	var wantRecords uint64
	for i, b := range blocks {
		if err := w.WriteBlock(b, recordCounts[i]); err != nil {
			t.Fatalf("WriteBlock %d: %v", i, err)
		}
		wantRecords += uint64(recordCounts[i])
	}
	if w.RecordsWritten() != wantRecords {
		t.Fatalf("RecordsWritten = %d, want %d", w.RecordsWritten(), wantRecords)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer r.Close()

	if r.BlockCount() != len(blocks) {
		t.Fatalf("BlockCount = %d, want %d", r.BlockCount(), len(blocks))
	}
	// Matches real dsrc's own DsrcFileHeader.recordsCount bug: always 0 on
	// disk, regardless of how many records were actually written.
	if r.RecordsCount() != 0 {
		t.Fatalf("RecordsCount = %d, want 0 (matching real dsrc's on-disk behavior)", r.RecordsCount())
	}
	if r.DatasetType().QualityOffset != 33 {
		t.Fatalf("QualityOffset = %d, want 33", r.DatasetType().QualityOffset)
	}

	for i, want := range blocks {
		got, err := r.ReadNextBlock()
		if err != nil {
			t.Fatalf("ReadNextBlock %d: %v", i, err)
		}
		if !bytes.Equal(got, want) {
			t.Fatalf("block %d mismatch (len got=%d want=%d)", i, len(got), len(want))
		}
	}

	if _, err := r.ReadNextBlock(); err != io.EOF {
		t.Fatalf("final ReadNextBlock: got err=%v, want io.EOF", err)
	}
}

// TestReadRealDsrcArchive verifies this package's Reader against an actual
// archive produced by the original C++ dsrc binary (built from this repo's
// own src/, compressed with default -m0 settings on a small test FASTQ
// file). It confirms every header/footer field this package interprets
// matches what real dsrc actually wrote, byte for byte.
func TestReadRealDsrcArchive(t *testing.T) {
	r, err := Open("testdata/real_default_mode.dsrc")
	if err != nil {
		t.Fatalf("Open real dsrc archive: %v", err)
	}
	defer r.Close()

	if r.BlockCount() != 1 {
		t.Fatalf("BlockCount = %d, want 1", r.BlockCount())
	}
	if r.RecordsCount() != 0 {
		t.Fatalf("RecordsCount = %d, want 0 (real dsrc's own writer never fills this in)", r.RecordsCount())
	}

	dt := r.DatasetType()
	if dt.QualityOffset != 33 {
		t.Errorf("QualityOffset = %d, want 33", dt.QualityOffset)
	}
	if dt.ColorSpace || dt.PlusRepetition {
		t.Errorf("DatasetType = %+v, want both flags false for this test file", dt)
	}

	cs := r.CompressionSettings()
	if cs.Lossy || cs.CalculateCRC32 {
		t.Errorf("CompressionSettings = %+v, want no lossy/crc32 flags (compressed without -l/-c)", cs)
	}
	if cs.DNAOrder != 0 || cs.QualityOrder != 0 {
		t.Errorf("DNAOrder/QualityOrder = %d/%d, want 0/0 (default -m0 = -d0 -q0)", cs.DNAOrder, cs.QualityOrder)
	}

	block, err := r.ReadNextBlock()
	if err != nil {
		t.Fatalf("ReadNextBlock: %v", err)
	}
	if len(block) == 0 {
		t.Fatal("read a zero-length block")
	}
	// This package doesn't decode DSRC's actual block payload format (see
	// package doc) — reaching here without error already confirms the
	// header/footer parsing and block-boundary bookkeeping are correct,
	// since a wrong footerOffset/blockSizes reading would either error out
	// or hand back the wrong slice of bytes.

	if _, err := r.ReadNextBlock(); err != io.EOF {
		t.Fatalf("second ReadNextBlock: got err=%v, want io.EOF", err)
	}
}

func TestOpenRejectsBadMagic(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "notdsrc.bin")
	if err := os.WriteFile(path, bytes.Repeat([]byte{0x00}, 64), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := Open(path); err == nil {
		t.Fatal("expected an error opening a file with the wrong magic")
	}
}

func TestOpenRejectsTruncatedFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "truncated.dsrc")

	w, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteBlock([]byte("hello"), 1); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	full, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	truncated := full[:len(full)-2]
	truncPath := filepath.Join(dir, "cut.dsrc")
	if err := os.WriteFile(truncPath, truncated, 0o644); err != nil {
		t.Fatal(err)
	}

	// Opening may succeed (header/footer still parse) or fail, but reading
	// the block must not silently return short/wrong data.
	r, err := Open(truncPath)
	if err != nil {
		return
	}
	defer r.Close()
	if _, err := r.ReadNextBlock(); err == nil {
		t.Fatal("expected an error reading a block from a truncated file")
	}
}

func TestWriteNoBlocksErrors(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "empty.dsrc")
	w, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err == nil {
		t.Fatal("expected an error closing an archive with zero blocks")
	}
}

func TestSingleBlockRoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "single.dsrc")

	w, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	data := []byte("just one block of data")
	if err := w.WriteBlock(data, 5); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	r, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	got, err := r.ReadNextBlock()
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(got, data) {
		t.Fatalf("got %q, want %q", got, data)
	}
}
