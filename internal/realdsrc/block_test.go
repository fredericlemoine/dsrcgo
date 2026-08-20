package realdsrc

import (
	"bytes"
	"os"
	"testing"
)

// TestBlockMatchesRealArchive verifies EncodeBlock/DecodeBlock end to end
// against a real dsrc archive: 100 records with varying read lengths
// (20/30/40/50 bases), which exercises the variable-length quality-length
// interleaving in tags (see tag.go's EncodeTags doc) and the resulting
// larger chunk-metadata section (StoreMetaData writes an extra minQuaLength
// word when FLAG_VARIABLE_LENGTH is set).
func TestBlockMatchesRealArchive(t *testing.T) {
	wantBlock := realBlock(t, "testdata/qvarlen.dsrc")

	chunkBytes, err := os.ReadFile("testdata/qvarlen.fastq")
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	// The archive's block payload doesn't include the file's own trailing
	// newline (see go/internal/fastq's ChunkReader, which trims it the
	// same way real dsrc's own chunk reader does).
	chunkBytes = bytes.TrimSuffix(chunkBytes, []byte("\n"))

	got, err := EncodeBlock(chunkBytes, 33)
	if err != nil {
		t.Fatalf("EncodeBlock: %v", err)
	}

	if !bytes.Equal(got, wantBlock) {
		t.Fatalf("block bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantBlock), got, wantBlock)
	}

	// DecodeBlock reconstructs a trailing newline after every record
	// (including the last), matching what a real FASTQ file looks like —
	// so compare against the untrimmed original, not the EncodeBlock input.
	original, err := os.ReadFile("testdata/qvarlen.fastq")
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}

	decoded, err := DecodeBlock(wantBlock, 33, false)
	if err != nil {
		t.Fatalf("DecodeBlock: %v", err)
	}
	if !bytes.Equal(decoded, original) {
		t.Fatalf("decoded chunk mismatch\n got  %q\n want %q", decoded, original)
	}
}

// TestBlockRoundTripPureACGTFixedLength exercises EncodeBlock/DecodeBlock
// on the earlier mini_2records fixture (pure ACGT, fixed length, B2 DNA
// scheme, Plain quality scheme) as a full end-to-end sanity check.
func TestBlockRoundTripPureACGTFixedLength(t *testing.T) {
	wantBlock := realBlock(t, "testdata/mini_2records.dsrc")

	// mini.fastq's raw bytes aren't saved as a fixture; reconstruct the
	// exact chunk text that produced mini_2records.dsrc directly.
	chunk := []byte("@R1\nACGTACGT\n+\nIIIIIIII\n@R2\nGGCCTTAA\n+\nIIIIIIII")

	got, err := EncodeBlock(chunk, 33)
	if err != nil {
		t.Fatalf("EncodeBlock: %v", err)
	}
	if !bytes.Equal(got, wantBlock) {
		t.Fatalf("block bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantBlock), got, wantBlock)
	}

	decoded, err := DecodeBlock(wantBlock, 33, false)
	if err != nil {
		t.Fatalf("DecodeBlock: %v", err)
	}
	// DecodeBlock adds back the trailing newline EncodeBlock's input
	// omitted (see the other test in this file for why).
	if !bytes.Equal(decoded, append(chunk, '\n')) {
		t.Fatalf("decoded chunk mismatch\n got  %q\n want %q", decoded, append(chunk, '\n'))
	}
}
