package dna

import "testing"

func TestSymbolTableRoundTrip(t *testing.T) {
	seq := []byte("ACGTNRWSKMDVHBYXU.-ACGT")
	indices, ok := EncodeSequence(seq)
	if !ok {
		t.Fatal("EncodeSequence: unexpected invalid byte")
	}
	got := DecodeSequence(indices)
	if string(got) != string(seq) {
		t.Fatalf("got %q, want %q", got, seq)
	}
}

func TestSymbolTableKnownIndices(t *testing.T) {
	cases := map[byte]byte{'A': 0, 'G': 1, 'C': 2, 'T': 3, 'N': 4, '-': 18}
	for b, want := range cases {
		got, ok := ToIndex(b)
		if !ok || got != want {
			t.Errorf("ToIndex(%q) = (%d, %v), want (%d, true)", b, got, ok, want)
		}
	}
}

func TestSymbolTableRejectsUnknownByte(t *testing.T) {
	if _, ok := ToIndex('z'); ok {
		t.Error("ToIndex('z') = ok, want rejected (lowercase not in IUPAC table)")
	}
	if _, ok := EncodeSequence([]byte("ACGTz")); ok {
		t.Error("EncodeSequence with invalid byte should fail")
	}
}
