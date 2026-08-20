package fastq

import "testing"

const sampleFastq = "@SEQ_ID_1 desc\n" +
	"GATTACAGATTACA\n" +
	"+SEQ_ID_1 desc\n" +
	"IIIIIIIIIIIIII\n" +
	"@SEQ_ID_2 desc\n" +
	"ACGTACGTACGTAC\n" +
	"+SEQ_ID_2 desc\n" +
	"HHHHHHHHHHHHHH\n"

func TestParseChunk(t *testing.T) {
	records, sizes, err := ParseChunk([]byte(sampleFastq))
	if err != nil {
		t.Fatalf("ParseChunk: %v", err)
	}
	if len(records) != 2 {
		t.Fatalf("got %d records, want 2", len(records))
	}

	r0 := records[0]
	if string(r0.Title) != "@SEQ_ID_1 desc" {
		t.Errorf("record 0 title = %q", r0.Title)
	}
	if string(r0.Sequence) != "GATTACAGATTACA" {
		t.Errorf("record 0 sequence = %q", r0.Sequence)
	}
	if string(r0.Quality) != "IIIIIIIIIIIIII" {
		t.Errorf("record 0 quality = %q", r0.Quality)
	}

	r1 := records[1]
	if string(r1.Title) != "@SEQ_ID_2 desc" {
		t.Errorf("record 1 title = %q", r1.Title)
	}

	wantTag := uint64(len(r0.Title) + len(r1.Title))
	if sizes.Tag != wantTag {
		t.Errorf("Tag size = %d, want %d", sizes.Tag, wantTag)
	}
	wantDna := uint64(len(r0.Sequence) + len(r1.Sequence))
	if sizes.Dna != wantDna {
		t.Errorf("Dna size = %d, want %d", sizes.Dna, wantDna)
	}
}

func TestParseChunkCRLF(t *testing.T) {
	data := "@id\r\nACGT\r\n+\r\nIIII\r\n@id2\r\nACGT\r\n+\r\nIIII\r\n"
	records, _, err := ParseChunk([]byte(data))
	if err != nil {
		t.Fatalf("ParseChunk: %v", err)
	}
	if len(records) != 2 {
		t.Fatalf("got %d records, want 2", len(records))
	}
	if string(records[0].Sequence) != "ACGT" {
		t.Errorf("sequence = %q", records[0].Sequence)
	}
}

func TestParseChunkEmpty(t *testing.T) {
	_, _, err := ParseChunk(nil)
	if err != ErrNoRecords {
		t.Fatalf("got err %v, want ErrNoRecords", err)
	}
}

func TestParseChunkTruncatedTrailingRecord(t *testing.T) {
	// A chunk cut mid-record (no trailing quality line) should still yield
	// the complete leading records, matching upstream's "stop at first
	// invalid record" behavior.
	data := sampleFastq + "@SEQ_ID_3 desc\nACGT\n+\n"
	records, _, err := ParseChunk([]byte(data))
	if err != nil {
		t.Fatalf("ParseChunk: %v", err)
	}
	if len(records) != 2 {
		t.Fatalf("got %d records, want 2 (truncated 3rd record should be dropped)", len(records))
	}
}

func TestAnalyzeQualityOffset(t *testing.T) {
	ds, ok := Analyze([]byte(sampleFastq), true)
	if !ok {
		t.Fatal("Analyze returned ok=false")
	}
	if ds.QualityOffset != 33 {
		t.Errorf("QualityOffset = %d, want 33", ds.QualityOffset)
	}
	if ds.ColorSpace {
		t.Error("ColorSpace = true, want false")
	}
	if !ds.PlusRepetition {
		t.Error("PlusRepetition = false, want true (record 1's + line repeats the tag)")
	}
}

func TestAnalyzeTooFewRecords(t *testing.T) {
	data := "@only_one\nACGT\n+\nIIII\n"
	_, ok := Analyze([]byte(data), false)
	if ok {
		t.Error("Analyze returned ok=true for a single-record chunk, want false")
	}
}
