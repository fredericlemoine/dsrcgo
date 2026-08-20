package huffman

import (
	"math/rand"
	"testing"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
)

// buildAndSerialize inserts freqs into a fresh tree, completes it, and
// returns both the codes (for encoding) and a freshly reloaded decode-side
// tree built purely from the serialized bitstream, mirroring how DSRC's
// encoder and decoder each end up with their own tree instance.
func buildAndSerialize(t *testing.T, freqs []uint32) (codes []Code, loaded *Tree) {
	t.Helper()

	enc := NewTree()
	enc.Restart(len(freqs))
	for _, f := range freqs {
		if !enc.Insert(f) {
			t.Fatal("Insert failed unexpectedly")
		}
	}
	codes = enc.Complete()
	if len(codes) != len(freqs) {
		t.Fatalf("Complete returned %d codes, want %d", len(codes), len(freqs))
	}

	w := bitio.NewWriter()
	enc.StoreTree(w)

	loaded = LoadTree(bitio.NewReader(w.Bytes()))
	return codes, loaded
}

func encodeSymbols(codes []Code, symbols []int) []byte {
	w := bitio.NewWriter()
	for _, s := range symbols {
		c := codes[s]
		w.PutBits(c.Code, uint(c.Len))
	}
	return w.Bytes()
}

func decodeSymbols(t *testing.T, tree *Tree, data []byte, n int) []int {
	t.Helper()
	r := bitio.NewReader(data)
	out := make([]int, n)
	for i := range out {
		out[i] = int(tree.DecodeSymbol(r))
	}
	return out
}

func TestTreeRoundTripSkewedFrequencies(t *testing.T) {
	freqs := []uint32{1000, 500, 200, 50, 10, 5, 1, 1}
	codes, loaded := buildAndSerialize(t, freqs)

	rng := rand.New(rand.NewSource(1))
	symbols := make([]int, 20000)
	for i := range symbols {
		// Sample roughly proportional to freqs.
		x := rng.Intn(1767) // sum(freqs)
		acc := 0
		for s, f := range freqs {
			acc += int(f)
			if x < acc {
				symbols[i] = s
				break
			}
		}
	}

	data := encodeSymbols(codes, symbols)
	got := decodeSymbols(t, loaded, data, len(symbols))

	for i := range symbols {
		if got[i] != symbols[i] {
			t.Fatalf("symbol %d: got %d, want %d", i, got[i], symbols[i])
		}
	}

	// A skewed distribution should compress well below the naive 3-bit
	// (8 symbols) packing.
	naiveBits := len(symbols) * 3
	gotBits := len(data) * 8
	if gotBits >= naiveBits {
		t.Errorf("didn't beat naive packing: %d bits vs naive %d bits", gotBits, naiveBits)
	}
	t.Logf("%d symbols -> %d bytes (%.3f bits/symbol, naive=3.0)", len(symbols), len(data), float64(gotBits)/float64(len(symbols)))
}

func TestTreeRoundTripUniformFrequencies(t *testing.T) {
	freqs := []uint32{10, 10, 10, 10, 10, 10, 10, 10, 10, 10}
	codes, loaded := buildAndSerialize(t, freqs)

	rng := rand.New(rand.NewSource(2))
	symbols := make([]int, 5000)
	for i := range symbols {
		symbols[i] = rng.Intn(len(freqs))
	}

	data := encodeSymbols(codes, symbols)
	got := decodeSymbols(t, loaded, data, len(symbols))
	for i := range symbols {
		if got[i] != symbols[i] {
			t.Fatalf("symbol %d: got %d, want %d", i, got[i], symbols[i])
		}
	}
}

func TestTreeSingleSymbol(t *testing.T) {
	// Only one distinct symbol ever appears — the n_symbols<2 special case.
	codes, loaded := buildAndSerialize(t, []uint32{42})

	data := encodeSymbols(codes, []int{0, 0, 0, 0})
	got := decodeSymbols(t, loaded, data, 4)
	for i, s := range got {
		if s != 0 {
			t.Fatalf("symbol %d: got %d, want 0", i, s)
		}
	}
}

func TestTreeWithZeroFrequencySymbols(t *testing.T) {
	// Symbols present in the alphabet but never observed at this position —
	// the common case for per-position quality trees. They must get a
	// valid (if unused) code and never be looked up.
	freqs := []uint32{100, 0, 0, 50, 0}
	codes, loaded := buildAndSerialize(t, freqs)

	symbols := []int{0, 3, 0, 3, 0, 0, 3}
	data := encodeSymbols(codes, symbols)
	got := decodeSymbols(t, loaded, data, len(symbols))
	for i := range symbols {
		if got[i] != symbols[i] {
			t.Fatalf("symbol %d: got %d, want %d", i, got[i], symbols[i])
		}
	}
}

func TestTreeManyPositionsIndependent(t *testing.T) {
	// Mirrors real usage: one independent Tree per read position, each
	// with its own frequency distribution.
	rng := rand.New(rand.NewSource(4))
	const positions = 100
	const alphabet = 8

	allCodes := make([][]Code, positions)
	allTrees := make([]*Tree, positions)
	freqsPerPos := make([][]uint32, positions)

	for p := 0; p < positions; p++ {
		freqs := make([]uint32, alphabet)
		for i := range freqs {
			freqs[i] = uint32(rng.Intn(1000))
		}
		freqsPerPos[p] = freqs
		allCodes[p], allTrees[p] = buildAndSerialize(t, freqs)
	}

	for p := 0; p < positions; p++ {
		symbols := make([]int, 200)
		for i := range symbols {
			symbols[i] = rng.Intn(alphabet)
		}
		data := encodeSymbols(allCodes[p], symbols)
		got := decodeSymbols(t, allTrees[p], data, len(symbols))
		for i := range symbols {
			if got[i] != symbols[i] {
				t.Fatalf("position %d symbol %d: got %d, want %d", p, i, got[i], symbols[i])
			}
		}
	}
}
