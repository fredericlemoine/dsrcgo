package bitio

import (
	"math/rand"
	"testing"
)

func TestPutGetBitsRoundTrip(t *testing.T) {
	rng := rand.New(rand.NewSource(1))

	type entry struct {
		word uint32
		n    uint
	}
	var entries []entry
	for i := 0; i < 5000; i++ {
		n := uint(1 + rng.Intn(32))
		var word uint32
		if n == 32 {
			word = rng.Uint32()
		} else {
			word = rng.Uint32() & ((1 << n) - 1)
		}
		entries = append(entries, entry{word, n})
	}

	w := NewWriter()
	for _, e := range entries {
		w.PutBits(e.word, e.n)
	}
	data := w.Bytes()

	r := NewReader(data)
	for i, e := range entries {
		got := r.GetBits(e.n)
		if got != e.word {
			t.Fatalf("entry %d: got %d, want %d (n=%d)", i, got, e.word, e.n)
		}
	}
}

func TestPutBitAndByteAndWord(t *testing.T) {
	w := NewWriter()
	w.PutBit(1)
	w.PutBit(0)
	w.PutByte(0xAB)
	w.PutWord(0xDEADBEEF)
	w.PutBit(1)

	r := NewReader(w.Bytes())
	if r.GetBit() != 1 {
		t.Fatal("bit 0 mismatch")
	}
	if r.GetBit() != 0 {
		t.Fatal("bit 1 mismatch")
	}
	if b := r.GetByte(); b != 0xAB {
		t.Fatalf("byte mismatch: got %#x", b)
	}
	if w32 := r.GetWord(); w32 != 0xDEADBEEF {
		t.Fatalf("word mismatch: got %#x", w32)
	}
	if r.GetBit() != 1 {
		t.Fatal("trailing bit mismatch")
	}
}

func TestZeroWidthIsNoOp(t *testing.T) {
	w := NewWriter()
	w.PutBits(0xFF, 0) // must not corrupt subsequent output
	w.PutByte(0x42)
	r := NewReader(w.Bytes())
	if r.GetBits(0) != 0 {
		t.Fatal("GetBits(0) should return 0")
	}
	if b := r.GetByte(); b != 0x42 {
		t.Fatalf("got %#x, want 0x42", b)
	}
}
