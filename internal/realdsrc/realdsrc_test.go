package realdsrc

import (
	"bytes"
	"math/rand"
	"testing"

	"github.com/fredericlemoine/dsrcgo/internal/archive"
	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/dna"
)

// realBlock reads the single block out of a real dsrc archive fixture.
func realBlock(t *testing.T, path string) []byte {
	t.Helper()
	r, err := archive.Open(path)
	if err != nil {
		t.Fatalf("archive.Open(%s): %v", path, err)
	}
	defer r.Close()
	block, err := r.ReadNextBlock()
	if err != nil {
		t.Fatalf("ReadNextBlock: %v", err)
	}
	return block
}

// indexStreams converts ASCII sequences to raw dna-table index streams —
// EncodeDNA/DecodeDNA work on the preprocessed representation (see
// preprocess.go), so tests exercising DNA in isolation (no smuggling)
// build that representation directly rather than through PreprocessForward.
func indexStreams(t *testing.T, sequences [][]byte) [][]byte {
	t.Helper()
	out := make([][]byte, len(sequences))
	for i, s := range sequences {
		idx, ok := dna.EncodeSequence(s)
		if !ok {
			t.Fatalf("sequence %d contains a non-IUPAC byte", i)
		}
		out[i] = idx
	}
	return out
}

func asciiStreams(indexed [][]byte) [][]byte {
	out := make([][]byte, len(indexed))
	for i, idx := range indexed {
		out[i] = dna.DecodeSequence(idx)
	}
	return out
}

// TestChunkHeaderMatchesRealArchive verifies WriteChunkHeader/ReadChunkHeader
// against the first 16 bytes of a real dsrc block: the mini_2records fixture
// is 2 records of length 8 with no variable-length flag, so
// recordsCount=2, maxQuaLength=8, flags=0, chunkSize=47 — values confirmed
// by manually decoding the real dsrc binary's own DEBUG-build output (which
// embeds section-boundary position markers) before writing this test.
func TestChunkHeaderMatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/mini_2records.dsrc")
	wantMeta := block[:16]

	h := ChunkHeader{RecordsCount: 2, MaxQuaLength: 8, MinQuaLength: 8, Flags: 0, ChunkSize: 47}
	w := bitio.NewWriter()
	WriteChunkHeader(w, h)
	got := w.Bytes()

	if !bytes.Equal(got, wantMeta) {
		t.Fatalf("metadata bytes mismatch\n got  %x\n want %x", got, wantMeta)
	}

	gotH := ReadChunkHeader(bitio.NewReader(wantMeta))
	if gotH != h {
		t.Fatalf("ReadChunkHeader = %+v, want %+v", gotH, h)
	}
}

// TestDNAB2MatchesRealArchive verifies EncodeDNA/DecodeDNA against the last
// 5 bytes of the same fixture: scheme byte (0x00 = B2) followed by 2
// records × 8 bases × 2 bits = 4 bytes of packed ACGT.
func TestDNAB2MatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/mini_2records.dsrc")
	wantDNA := block[len(block)-5:]

	sequences := [][]byte{[]byte("ACGTACGT"), []byte("GGCCTTAA")}
	streams := indexStreams(t, sequences)

	w := bitio.NewWriter()
	if err := EncodeDNA(w, streams); err != nil {
		t.Fatalf("EncodeDNA: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantDNA) {
		t.Fatalf("DNA bytes mismatch\n got  %x\n want %x", got, wantDNA)
	}

	decoded, err := DecodeDNA(bitio.NewReader(wantDNA), []int{8, 8})
	if err != nil {
		t.Fatalf("DecodeDNA: %v", err)
	}
	decodedASCII := asciiStreams(decoded)
	for i, want := range sequences {
		if !bytes.Equal(decodedASCII[i], want) {
			t.Fatalf("sequence %d: got %q, want %q", i, decodedASCII[i], want)
		}
	}
}

func TestDNAB2RoundTripRandomACGT(t *testing.T) {
	sequences := [][]byte{
		[]byte("AAAA"),
		[]byte("TTTTTTTT"),
		[]byte("ACGTACGTACGT"),
		[]byte("GGGGCCCC"),
	}
	streams := indexStreams(t, sequences)
	lengths := make([]int, len(streams))
	for i, s := range streams {
		lengths[i] = len(s)
	}

	w := bitio.NewWriter()
	if err := EncodeDNA(w, streams); err != nil {
		t.Fatalf("EncodeDNA: %v", err)
	}

	decoded, err := DecodeDNA(bitio.NewReader(w.Bytes()), lengths)
	if err != nil {
		t.Fatalf("DecodeDNA: %v", err)
	}
	decodedASCII := asciiStreams(decoded)
	for i, want := range sequences {
		if !bytes.Equal(decodedASCII[i], want) {
			t.Fatalf("sequence %d: got %q, want %q", i, decodedASCII[i], want)
		}
	}
}

func TestDNAB2RejectsCorruptingCase(t *testing.T) {
	// Only A, T, N appear (3 distinct symbols, <=4): real dsrc would still
	// pick SchemeB2 here since it only checks the count, then silently
	// truncate N's fixed index (4) to 0 via Put2Bits's 2-bit mask,
	// corrupting the data. This package refuses instead — see dna.go's
	// package doc.
	streams := indexStreams(t, [][]byte{[]byte("AATTNN")})
	w := bitio.NewWriter()
	err := EncodeDNA(w, streams)
	if err == nil {
		t.Fatal("expected an error for a <=4-symbol block that isn't pure ACGT")
	}
}

// TestDNAHuffmanMatchesRealArchive verifies EncodeDNA/DecodeDNA's Huffman
// path against a real dsrc archive: 5 records, each 10 bases of a single
// repeated symbol (A/G/C/T/N respectively, with one extra 'A' swapped into
// the N record), giving frequencies A=11,G=10,C=10,T=10,N=9 — skewed
// enough to produce a non-trivial, non-balanced tree. The expected tree
// shape (and thus these exact bytes) was hand-verified by decoding real
// dsrc's own DEBUG-build output bit by bit before this test was written.
func TestDNAHuffmanMatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/dna5_huffman.dsrc")
	wantDNA := block[len(block)-35:]

	sequences := [][]byte{
		[]byte("AAAAAAAAAA"),
		[]byte("GGGGGGGGGG"),
		[]byte("CCCCCCCCCC"),
		[]byte("TTTTTTTTTT"),
		[]byte("NNNNNNNNNA"),
	}
	streams := indexStreams(t, sequences)

	w := bitio.NewWriter()
	if err := EncodeDNA(w, streams); err != nil {
		t.Fatalf("EncodeDNA: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantDNA) {
		t.Fatalf("DNA bytes mismatch\n got  %x\n want %x", got, wantDNA)
	}

	lengths := make([]int, len(sequences))
	for i, s := range sequences {
		lengths[i] = len(s)
	}
	decoded, err := DecodeDNA(bitio.NewReader(wantDNA), lengths)
	if err != nil {
		t.Fatalf("DecodeDNA: %v", err)
	}
	decodedASCII := asciiStreams(decoded)
	for i, want := range sequences {
		if !bytes.Equal(decodedASCII[i], want) {
			t.Fatalf("sequence %d: got %q, want %q", i, decodedASCII[i], want)
		}
	}
}

func TestDNAHuffmanRoundTripLargerAlphabet(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	bases := []byte("AGCTNRWSKM") // 10 distinct symbols, forces the Huffman path
	sequences := make([][]byte, 200)
	for i := range sequences {
		n := 20 + rng.Intn(80)
		seq := make([]byte, n)
		for j := range seq {
			seq[j] = bases[rng.Intn(len(bases))]
		}
		sequences[i] = seq
	}
	streams := indexStreams(t, sequences)
	lengths := make([]int, len(streams))
	for i, s := range streams {
		lengths[i] = len(s)
	}

	w := bitio.NewWriter()
	if err := EncodeDNA(w, streams); err != nil {
		t.Fatalf("EncodeDNA: %v", err)
	}

	decoded, err := DecodeDNA(bitio.NewReader(w.Bytes()), lengths)
	if err != nil {
		t.Fatalf("DecodeDNA: %v", err)
	}
	decodedASCII := asciiStreams(decoded)
	for i, want := range sequences {
		if !bytes.Equal(decodedASCII[i], want) {
			t.Fatalf("sequence %d: got %q, want %q", i, decodedASCII[i], want)
		}
	}
}
