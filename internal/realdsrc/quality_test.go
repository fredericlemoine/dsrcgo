package realdsrc

import (
	"bytes"
	"math/rand"
	"os"
	"testing"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/fastq"
)

// qualityStreamsFromFile parses a FASTQ file and returns each record's
// quality string converted to 0-based numeric scores (ASCII - offset) —
// the same representation PreprocessForward's .Quality field holds when no
// ambiguous-base smuggling occurs (pure ACGT input, as in the fixture
// below), so it can be fed to EncodeQuality directly.
func qualityStreamsFromFile(t *testing.T, path string, offset byte) [][]byte {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile(%s): %v", path, err)
	}
	records, _, err := fastq.ParseChunk(data)
	if err != nil {
		t.Fatalf("ParseChunk: %v", err)
	}
	out := make([][]byte, len(records))
	for i, r := range records {
		q := make([]byte, len(r.Quality))
		for j, b := range r.Quality {
			q[j] = b - offset
		}
		out[i] = q
	}
	return out
}

// TestQualityPlainMatchesRealArchive verifies EncodeQuality/DecodeQuality's
// Plain scheme against a real dsrc archive: 30 records of 20 bases each
// with realistic, non-repetitive quality scores (a decaying random walk),
// chosen specifically to avoid triggering the RLE scheme (needs long runs
// of a repeated value) or the Truncated scheme (needs a long low-quality
// tail) — real dsrc's own scheme-selection byte in this archive confirms
// it picked Plain (verified via the debug build before writing this test).
// Section boundaries (metadata=16B, tags=98B, quality=807B, dna=151B) were
// also determined from the debug build's embedded position markers.
func TestQualityPlainMatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/qplain.dsrc")
	const metaLen, tagsLen, qualityLen = 16, 98, 807
	wantQuality := block[metaLen+tagsLen : metaLen+tagsLen+qualityLen]

	streams := qualityStreamsFromFile(t, "testdata/qplain.fastq", 33)
	if len(streams) != 30 {
		t.Fatalf("got %d records, want 30", len(streams))
	}

	st := computeQualityStats(streams)
	if scheme := selectQualityScheme(st); scheme != qualitySchemePlain {
		t.Fatalf("selectQualityScheme = %d, want Plain (0) — fixture no longer matches real dsrc's own scheme choice", scheme)
	}

	w := bitio.NewWriter()
	if err := EncodeQuality(w, streams); err != nil {
		t.Fatalf("EncodeQuality: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantQuality) {
		t.Fatalf("quality bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantQuality), got, wantQuality)
	}

	lengths := make([]int, len(streams))
	for i, s := range streams {
		lengths[i] = len(s)
	}
	decoded, err := DecodeQuality(bitio.NewReader(wantQuality), lengths)
	if err != nil {
		t.Fatalf("DecodeQuality: %v", err)
	}
	for i, want := range streams {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("record %d: got %v, want %v", i, decoded[i], want)
		}
	}
}

// TestQualityTruncatedMatchesRealArchive verifies EncodeQuality/
// DecodeQuality's Truncated scheme against a real dsrc archive: 30 records
// of 20 bases, each with 12 positions of varied quality followed by a
// trailing run of 8 Phred-2 ("poor") calls — the trailing run triggers
// real dsrc's rawLength/thLength > 1.10 truncation heuristic (confirmed
// via the debug build's scheme-selection byte, 1 = Truncated). Section
// boundaries (metadata=16B, tags=95B, quality=660B, dna=151B) were also
// determined from the debug build's embedded position markers.
func TestQualityTruncatedMatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/qtrunc.dsrc")
	const metaLen, tagsLen, qualityLen = 16, 95, 660
	wantQuality := block[metaLen+tagsLen : metaLen+tagsLen+qualityLen]

	streams := qualityStreamsFromFile(t, "testdata/qtrunc.fastq", 33)
	if len(streams) != 30 {
		t.Fatalf("got %d records, want 30", len(streams))
	}

	st := computeQualityStats(streams)
	if scheme := selectQualityScheme(st); scheme != qualitySchemeTruncated {
		t.Fatalf("selectQualityScheme = %d, want Truncated (1) — fixture no longer matches real dsrc's own scheme choice", scheme)
	}

	w := bitio.NewWriter()
	if err := EncodeQuality(w, streams); err != nil {
		t.Fatalf("EncodeQuality: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantQuality) {
		t.Fatalf("quality bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantQuality), got, wantQuality)
	}

	lengths := make([]int, len(streams))
	for i, s := range streams {
		lengths[i] = len(s)
	}
	decoded, err := DecodeQuality(bitio.NewReader(wantQuality), lengths)
	if err != nil {
		t.Fatalf("DecodeQuality: %v", err)
	}
	for i, want := range streams {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("record %d: got %v, want %v", i, decoded[i], want)
		}
	}
}

func TestQualityTruncatedRoundTripSynthetic(t *testing.T) {
	rng := rand.New(rand.NewSource(11))
	streams := make([][]byte, 40)
	for i := range streams {
		s := make([]byte, 25)
		cur := 30
		for j := 0; j < 15; j++ {
			cur += rng.Intn(5) - 2
			if cur < 10 {
				cur = 10
			}
			if cur > 40 {
				cur = 40
			}
			s[j] = byte(cur)
		}
		for j := 15; j < 25; j++ {
			s[j] = hashSymbolNormalQ
		}
		streams[i] = s
	}

	st := computeQualityStats(streams)
	if scheme := selectQualityScheme(st); scheme != qualitySchemeTruncated {
		t.Fatalf("selectQualityScheme = %d, want Truncated (1)", scheme)
	}

	w := bitio.NewWriter()
	if err := EncodeQuality(w, streams); err != nil {
		t.Fatalf("EncodeQuality: %v", err)
	}

	lengths := make([]int, len(streams))
	for i, s := range streams {
		lengths[i] = len(s)
	}
	decoded, err := DecodeQuality(bitio.NewReader(w.Bytes()), lengths)
	if err != nil {
		t.Fatalf("DecodeQuality: %v", err)
	}
	for i, want := range streams {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("record %d: got %v, want %v", i, decoded[i], want)
		}
	}
}

func TestTruncatedLen(t *testing.T) {
	cases := []struct {
		q    []byte
		want int
	}{
		{[]byte{}, 0},
		{[]byte{5, 6, 7}, 3},
		{[]byte{5, 6, 2, 2, 2}, 2},
		{[]byte{2, 2, 2}, 1},
		{[]byte{2, 5, 2}, 2},
	}
	for _, c := range cases {
		if got := truncatedLen(c.q); got != c.want {
			t.Errorf("truncatedLen(%v) = %d, want %d", c.q, got, c.want)
		}
	}
}

func TestQualityPlainRoundTripSyntheticVariableLengths(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	lengths := []int{36, 100, 75, 151, 10}
	streams := make([][]byte, len(lengths))
	for i, l := range lengths {
		s := make([]byte, l)
		cur := 35
		for j := range s {
			cur += rng.Intn(5) - 2
			if cur < 2 {
				cur = 2
			}
			if cur > 40 {
				cur = 40
			}
			s[j] = byte(cur)
		}
		streams[i] = s
	}

	st := computeQualityStats(streams)
	scheme := selectQualityScheme(st)

	w := bitio.NewWriter()
	err := EncodeQuality(w, streams)
	if scheme != qualitySchemePlain {
		if err == nil {
			t.Fatalf("expected an error: this synthetic data selected scheme %d, which isn't implemented", scheme)
		}
		return
	}
	if err != nil {
		t.Fatalf("EncodeQuality: %v", err)
	}

	decoded, err := DecodeQuality(bitio.NewReader(w.Bytes()), lengths)
	if err != nil {
		t.Fatalf("DecodeQuality: %v", err)
	}
	for i, want := range streams {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("record %d: got %v, want %v", i, decoded[i], want)
		}
	}
}

func TestQualitySchemeSelectionRLE(t *testing.T) {
	// Long runs of a single repeated quality value strongly favor RLE.
	streams := make([][]byte, 50)
	for i := range streams {
		streams[i] = bytes.Repeat([]byte{40}, 100)
	}
	st := computeQualityStats(streams)
	if scheme := selectQualityScheme(st); scheme != qualitySchemeRLE {
		t.Fatalf("selectQualityScheme = %d, want RLE (2) for highly repetitive data", scheme)
	}
}

// TestQualityRLEMatchesRealArchive verifies EncodeQuality/DecodeQuality's
// RLE scheme (the context-conditional-Huffman path, qSymbolCount > 1)
// against a real dsrc archive: 20 records of 30 bases with quality built
// from long runs (5-15 each) of 4 distinct values, which drives
// thLength/rleLength past real dsrc's 1.25 RLE threshold while still
// having several distinct symbols to condition on. Real dsrc's own
// scheme-selection byte in this archive (2 = RLE) was confirmed via the
// debug build. Section boundaries (metadata=16B, tags=86B, quality=251B)
// were also determined from the debug build's embedded position markers.
func TestQualityRLEMatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/qrle.dsrc")
	const metaLen, tagsLen, qualityLen = 16, 86, 251
	wantQuality := block[metaLen+tagsLen : metaLen+tagsLen+qualityLen]

	streams := qualityStreamsFromFile(t, "testdata/qrle.fastq", 33)
	if len(streams) != 20 {
		t.Fatalf("got %d records, want 20", len(streams))
	}

	st := computeQualityStats(streams)
	if scheme := selectQualityScheme(st); scheme != qualitySchemeRLE {
		t.Fatalf("selectQualityScheme = %d, want RLE (2) — fixture no longer matches real dsrc's own scheme choice", scheme)
	}

	w := bitio.NewWriter()
	if err := EncodeQuality(w, streams); err != nil {
		t.Fatalf("EncodeQuality: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantQuality) {
		t.Fatalf("quality bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantQuality), got, wantQuality)
	}

	lengths := make([]int, len(streams))
	for i, s := range streams {
		lengths[i] = len(s)
	}
	decoded, err := DecodeQuality(bitio.NewReader(wantQuality), lengths)
	if err != nil {
		t.Fatalf("DecodeQuality: %v", err)
	}
	for i, want := range streams {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("record %d: got %v, want %v", i, decoded[i], want)
		}
	}
}

// TestQualityRLESingleSymbolMatchesRealArchive verifies the degenerate
// qSymbolCount==1 path: 20 records of 40 identical-quality bases (800
// total positions), forcing multiple runs via the 255-length cap without
// any symbol diversity to condition on. Real dsrc's scheme byte (2 = RLE)
// and section boundaries (metadata=16B, tags=86B, quality=70B) were
// confirmed via the debug build.
func TestQualityRLESingleSymbolMatchesRealArchive(t *testing.T) {
	block := realBlock(t, "testdata/qrle_single.dsrc")
	const metaLen, tagsLen, qualityLen = 16, 86, 70
	wantQuality := block[metaLen+tagsLen : metaLen+tagsLen+qualityLen]

	streams := qualityStreamsFromFile(t, "testdata/qrle_single.fastq", 33)
	if len(streams) != 20 {
		t.Fatalf("got %d records, want 20", len(streams))
	}

	st := computeQualityStats(streams)
	if scheme := selectQualityScheme(st); scheme != qualitySchemeRLE {
		t.Fatalf("selectQualityScheme = %d, want RLE (2)", scheme)
	}

	w := bitio.NewWriter()
	if err := EncodeQuality(w, streams); err != nil {
		t.Fatalf("EncodeQuality: %v", err)
	}
	got := w.Bytes()

	if !bytes.Equal(got, wantQuality) {
		t.Fatalf("quality bytes mismatch (len got=%d want=%d)\n got  %x\n want %x", len(got), len(wantQuality), got, wantQuality)
	}

	lengths := make([]int, len(streams))
	for i, s := range streams {
		lengths[i] = len(s)
	}
	decoded, err := DecodeQuality(bitio.NewReader(wantQuality), lengths)
	if err != nil {
		t.Fatalf("DecodeQuality: %v", err)
	}
	for i, want := range streams {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("record %d: got %v, want %v", i, decoded[i], want)
		}
	}
}

func TestQualityRLERoundTripSyntheticVariableLengths(t *testing.T) {
	rng := rand.New(rand.NewSource(13))
	symbols := []byte{5, 20, 35, 40}
	lengths := []int{36, 100, 75, 151}
	streams := make([][]byte, len(lengths))
	for i, l := range lengths {
		s := make([]byte, 0, l)
		for len(s) < l {
			sym := symbols[rng.Intn(len(symbols))]
			run := 5 + rng.Intn(20)
			for k := 0; k < run && len(s) < l; k++ {
				s = append(s, sym)
			}
		}
		streams[i] = s
	}

	st := computeQualityStats(streams)
	scheme := selectQualityScheme(st)

	w := bitio.NewWriter()
	if err := EncodeQuality(w, streams); err != nil {
		t.Fatalf("EncodeQuality (scheme %d): %v", scheme, err)
	}

	decoded, err := DecodeQuality(bitio.NewReader(w.Bytes()), lengths)
	if err != nil {
		t.Fatalf("DecodeQuality: %v", err)
	}
	for i, want := range streams {
		if !bytes.Equal(decoded[i], want) {
			t.Fatalf("record %d: got %v, want %v", i, decoded[i], want)
		}
	}
}
