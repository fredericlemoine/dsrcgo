package realdsrc

import (
	"bytes"
	"os"
	"testing"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/fastq"
)

func tagsFromFile(t *testing.T, path string) [][]byte {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile(%s): %v", path, err)
	}
	records, _, err := fastq.ParseChunk(data)
	if err != nil {
		t.Fatalf("ParseChunk: %v", err)
	}
	tags := make([][]byte, len(records))
	for i, r := range records {
		tags[i] = r.Title
	}
	return tags
}

// TestTagsMatchRealArchive verifies EncodeTags/DecodeTags against a real
// dsrc archive: 200 records with a 3-field tag
// "@SRR001471.<incrementing> HWI-EAS7:1:3:<x>:<y>/1" — the first numeric
// field is a pure incrementing counter (DeltaConst), x and y are random
// within a wide range (ValueVar, fixed-width). Section boundaries
// (metadata=16B, tags=804B) were determined from the debug build's
// embedded position markers.
func TestTagsMatchRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/qtags.dsrc")
	const metaLen, tagsLen = 16, 804
	wantTags := block[metaLen : metaLen+tagsLen]

	tags := tagsFromFile(t, "testdata/qtags.fastq")
	if len(tags) != 200 {
		t.Fatalf("got %d tags, want 200", len(tags))
	}

	w := bitio.NewWriter()
	if err := EncodeTags(w, tags, nil, 0, 0); err != nil {
		t.Fatalf("EncodeTags: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantTags) {
		t.Fatalf("tags bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantTags), got, wantTags)
	}

	decoded, _, err := DecodeTags(bitio.NewReader(wantTags), len(tags), 0, 0)
	if err != nil {
		t.Fatalf("DecodeTags: %v", err)
	}
	for i, want := range tags {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("tag %d: got %q, want %q", i, decoded[i], want)
		}
	}
}

// TestTagsTextFieldMatchesRealArchive verifies EncodeTags/DecodeTags
// against a real dsrc archive whose only non-constant field is free text
// (a 6-character barcode drawn from ACGTN, no numeric interpretation) —
// the hamming-mask + per-position Huffman path. Section boundaries
// (metadata=16B, tags=416B) were determined from the debug build's
// embedded position markers.
func TestTagsTextFieldMatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/qtagtext.dsrc")
	const metaLen, tagsLen = 16, 416
	wantTags := block[metaLen : metaLen+tagsLen]

	tags := tagsFromFile(t, "testdata/qtagtext.fastq")
	if len(tags) != 150 {
		t.Fatalf("got %d tags, want 150", len(tags))
	}

	w := bitio.NewWriter()
	if err := EncodeTags(w, tags, nil, 0, 0); err != nil {
		t.Fatalf("EncodeTags: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantTags) {
		t.Fatalf("tags bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantTags), got, wantTags)
	}

	decoded, _, err := DecodeTags(bitio.NewReader(wantTags), len(tags), 0, 0)
	if err != nil {
		t.Fatalf("DecodeTags: %v", err)
	}
	for i, want := range tags {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("tag %d: got %q, want %q", i, decoded[i], want)
		}
	}
}

// TestTagsValueVarHuffmanMatchesRealArchive verifies EncodeTags/DecodeTags
// against a real dsrc archive whose single numeric field has only 5
// distinct values spread across a wide range (100000..999999) with no two
// consecutive records repeating a value (avoiding the ValueRle scheme,
// which real dsrc picks instead when repeats cluster — confirmed
// empirically: an earlier version of this fixture that allowed consecutive
// repeats did select ValueRle). Real dsrc's own scheme byte here (1 =
// ValueVar) confirms the Huffman-optimized path is used, since only 5
// distinct values appear despite the wide range. Section boundaries
// (metadata=16B, tags=775B) were determined from the debug build's
// embedded position markers.
func TestTagsValueVarHuffmanMatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/qtagshuf.dsrc")
	const metaLen, tagsLen = 16, 775
	wantTags := block[metaLen : metaLen+tagsLen]

	tags := tagsFromFile(t, "testdata/qtagshuf.fastq")
	if len(tags) != 300 {
		t.Fatalf("got %d tags, want 300", len(tags))
	}

	w := bitio.NewWriter()
	if err := EncodeTags(w, tags, nil, 0, 0); err != nil {
		t.Fatalf("EncodeTags: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantTags) {
		t.Fatalf("tags bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantTags), got, wantTags)
	}

	decoded, _, err := DecodeTags(bitio.NewReader(wantTags), len(tags), 0, 0)
	if err != nil {
		t.Fatalf("DecodeTags: %v", err)
	}
	for i, want := range tags {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("tag %d: got %q, want %q", i, decoded[i], want)
		}
	}
}
