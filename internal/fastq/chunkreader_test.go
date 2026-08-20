package fastq

import (
	"bytes"
	"fmt"
	"io"
	"math/rand"
	"testing"
)

// buildFastq generates n synthetic FASTQ records with varying line lengths.
func buildFastq(n int, seed int64) []byte {
	rng := rand.New(rand.NewSource(seed))
	bases := "ACGT"
	var buf bytes.Buffer
	for i := 0; i < n; i++ {
		length := 20 + rng.Intn(80)
		seq := make([]byte, length)
		qual := make([]byte, length)
		for j := range seq {
			seq[j] = bases[rng.Intn(len(bases))]
			qual[j] = byte(33 + rng.Intn(40))
		}
		fmt.Fprintf(&buf, "@read%d some description here\n", i)
		buf.Write(seq)
		buf.WriteByte('\n')
		buf.WriteString("+\n")
		buf.Write(qual)
		buf.WriteByte('\n')
	}
	return buf.Bytes()
}

// TestChunkReaderRoundTrip verifies that splitting a synthetic FASTQ file
// into small chunks (forcing many boundary searches) and parsing each chunk
// reproduces exactly the records obtained by parsing the whole file at once.
func TestChunkReaderRoundTrip(t *testing.T) {
	data := buildFastq(5000, 42)

	wantRecords, _, err := ParseChunk(data)
	if err != nil {
		t.Fatalf("ParseChunk(whole file): %v", err)
	}

	// Small chunk size relative to record size forces the boundary-search
	// path (nextRecordBoundary) on nearly every chunk.
	const smallChunk = tailScanWindow + 512
	cr := NewChunkReader(bytes.NewReader(data), smallChunk)

	var gotRecords []Record
	chunkCount := 0
	for {
		chunk, err := cr.ReadNextChunk()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("ReadNextChunk: %v", err)
		}
		chunkCount++

		// Copy the chunk since ChunkReader reuses its internal buffer.
		owned := append([]byte(nil), chunk...)
		recs, _, err := ParseChunk(owned)
		if err != nil {
			t.Fatalf("ParseChunk(chunk %d): %v", chunkCount, err)
		}
		gotRecords = append(gotRecords, recs...)
	}

	if chunkCount < 2 {
		t.Fatalf("test is meaningless with only %d chunk(s); expected many", chunkCount)
	}

	if len(gotRecords) != len(wantRecords) {
		t.Fatalf("got %d records across %d chunks, want %d", len(gotRecords), chunkCount, len(wantRecords))
	}
	for i := range wantRecords {
		if !bytes.Equal(gotRecords[i].Title, wantRecords[i].Title) ||
			!bytes.Equal(gotRecords[i].Sequence, wantRecords[i].Sequence) ||
			!bytes.Equal(gotRecords[i].Quality, wantRecords[i].Quality) {
			t.Fatalf("record %d mismatch:\n got  title=%q seq=%q qual=%q\n want title=%q seq=%q qual=%q",
				i, gotRecords[i].Title, gotRecords[i].Sequence, gotRecords[i].Quality,
				wantRecords[i].Title, wantRecords[i].Sequence, wantRecords[i].Quality)
		}
	}
}

func TestChunkReaderSmallFile(t *testing.T) {
	data := buildFastq(3, 7)
	cr := NewChunkReader(bytes.NewReader(data), DefaultChunkSize)

	chunk, err := cr.ReadNextChunk()
	if err != nil {
		t.Fatalf("ReadNextChunk: %v", err)
	}
	recs, _, err := ParseChunk(chunk)
	if err != nil {
		t.Fatalf("ParseChunk: %v", err)
	}
	if len(recs) != 3 {
		t.Fatalf("got %d records, want 3", len(recs))
	}

	_, err = cr.ReadNextChunk()
	if err != io.EOF {
		t.Fatalf("second ReadNextChunk: got err=%v, want io.EOF", err)
	}
}
