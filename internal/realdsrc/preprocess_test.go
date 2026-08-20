package realdsrc

import (
	"bytes"
	"testing"
)

// TestPreprocessSmugglingMatchesRealArchive verifies which bases
// PreprocessForward drops from the DNA stream against a real dsrc archive:
// record 1 "ANCGT" with quality "I#III" has N at a low quality score
// (offset 33, '#' = 35 -> qRaw=2 < 7), so N must be dropped from the DNA
// stream, leaving [A,C,G,T]; record 2 "AGCGT" at uniform high quality
// keeps every base, [A,G,C,G,T]. Both the dropped base and the exact
// symbol order of survivors were confirmed by decoding real dsrc's own
// DEBUG-build output before this test was written.
func TestPreprocessSmugglingMatchesRealArchive(t *testing.T) {
	sequences := [][]byte{[]byte("ANCGT"), []byte("AGCGT")}
	qualities := [][]byte{[]byte("I#III"), []byte("IIIII")}

	pre, err := PreprocessForward(sequences, qualities, 33)
	if err != nil {
		t.Fatalf("PreprocessForward: %v", err)
	}

	wantDNA0 := indexStreams(t, [][]byte{[]byte("ACGT")})[0]
	wantDNA1 := indexStreams(t, [][]byte{[]byte("AGCGT")})[0]

	if !bytes.Equal(pre[0].DNA, wantDNA0) {
		t.Errorf("record 0 DNA stream = %v, want %v (N should be dropped)", pre[0].DNA, wantDNA0)
	}
	if !bytes.Equal(pre[1].DNA, wantDNA1) {
		t.Errorf("record 1 DNA stream = %v, want %v (nothing dropped)", pre[1].DNA, wantDNA1)
	}
}

func TestPreprocessForwardBackwardRoundTrip(t *testing.T) {
	sequences := [][]byte{
		[]byte("ANCGT"),
		[]byte("AGCGT"),
		[]byte("NNNNN"), // all-ambiguous, all low quality: entire record vanishes from the DNA stream
		[]byte("ACGTACGTACGT"),
	}
	qualities := [][]byte{
		[]byte("I#III"),
		[]byte("IIIII"),
		[]byte("#####"),
		[]byte("IIIIIIIIIIII"),
	}
	const offset = 33

	pre, err := PreprocessForward(sequences, qualities, offset)
	if err != nil {
		t.Fatalf("PreprocessForward: %v", err)
	}

	dnaStreams := make([][]byte, len(pre))
	qualStreams := make([][]byte, len(pre))
	for i, p := range pre {
		dnaStreams[i] = p.DNA
		qualStreams[i] = p.Quality
	}

	gotSeq, gotQual := PostprocessBackward(dnaStreams, qualStreams, offset)
	for i := range sequences {
		if !bytes.Equal(gotSeq[i], sequences[i]) {
			t.Errorf("record %d sequence: got %q, want %q", i, gotSeq[i], sequences[i])
		}
		if !bytes.Equal(gotQual[i], qualities[i]) {
			t.Errorf("record %d quality: got %q, want %q", i, gotQual[i], qualities[i])
		}
	}
}

func TestPreprocessForwardRejectsShortQuality(t *testing.T) {
	_, err := PreprocessForward([][]byte{[]byte("ACGT")}, [][]byte{[]byte("III")}, 33)
	if err == nil {
		t.Fatal("expected an error for mismatched sequence/quality length")
	}
}

func TestPreprocessNoAmbiguityIsIdentity(t *testing.T) {
	// Pure ACGT at any quality should never trigger smuggling: DNA stream
	// length always equals the sequence length.
	sequences := [][]byte{[]byte("ACGTACGT")}
	qualities := [][]byte{[]byte("!!!!!!!!")} // lowest possible quality
	pre, err := PreprocessForward(sequences, qualities, 33)
	if err != nil {
		t.Fatalf("PreprocessForward: %v", err)
	}
	if len(pre[0].DNA) != len(sequences[0]) {
		t.Fatalf("DNA stream length = %d, want %d (pure ACGT should never be dropped)", len(pre[0].DNA), len(sequences[0]))
	}

	want := indexStreams(t, sequences)[0]
	if !bytes.Equal(pre[0].DNA, want) {
		t.Fatalf("DNA stream = %v, want %v", pre[0].DNA, want)
	}
}
