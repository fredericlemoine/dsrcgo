// Package huffman ports DSRC's canonical Huffman coder (src/huffman.h,
// src/huffman.cpp) to Go. It's used by the default (order-0) quality
// compression scheme: one tree per read position, built from that
// position's symbol-frequency distribution.
//
// Two deliberate simplifications from upstream, both local to this package
// and invisible to callers:
//
//  1. Upstream keeps two different in-memory tree shapes — one built by
//     Complete() for code assignment (leaves are explicit nodes flagged by
//     left_child == -1) and a differently-shaped one rebuilt by LoadTree
//     for decoding (leaves are embedded as negative child pointers, no
//     separate node slot) — apparently to avoid allocating full leaf nodes
//     on the decode side. This port uses one tree shape for both, which
//     carries the same information and round-trips correctly; only the
//     internal bookkeeping differs.
//  2. DecodeFast's precomputed speedup table (src/huffman.cpp
//     ComputeSpeedupTree) is a pure performance optimization equivalent to
//     repeated single-bit Decode calls; this port only implements the
//     single-bit path.
//
// The single-symbol edge case (n_symbols forced from 1 to 2) also departs
// from upstream deliberately: upstream pads with a second heap slot that,
// for a tree Restart with capacity 1, reads past the end of its allocation
// — undefined behavior we do not reproduce. Here the padding symbol gets a
// fresh, in-bounds id instead; see Complete.
package huffman

import (
	"container/heap"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
)

// Code is a symbol's Huffman code: the low Len bits of Code, MSB first.
type Code struct {
	Code uint32
	Len  uint32
}

type node struct {
	left, right int32
}

type freqItem struct {
	symbol uint32
	freq   uint32
}

// freqHeap is a min-heap by frequency, ties broken toward the lower symbol
// id — matching the ordering HuffmanEncoder::Frequency::operator< induces
// via std::make_heap/pop_heap (see package doc for why: a max-heap under
// that inverted comparator is a min-heap by frequency).
type freqHeap []freqItem

func (h freqHeap) Len() int { return len(h) }
func (h freqHeap) Less(i, j int) bool {
	if h[i].freq != h[j].freq {
		return h[i].freq < h[j].freq
	}
	return h[i].symbol < h[j].symbol
}
func (h freqHeap) Swap(i, j int) { h[i], h[j] = h[j], h[i] }
func (h *freqHeap) Push(x any)   { *h = append(*h, x.(freqItem)) }
func (h *freqHeap) Pop() any {
	old := *h
	n := len(old)
	item := old[n-1]
	*h = old[:n-1]
	return item
}

// Tree is one Huffman coder instance: symbol frequencies go in via
// Insert/Complete to build encode codes, or a serialized shape goes in via
// LoadTree to enable decoding.
type Tree struct {
	capacity int
	freqs    []freqItem

	leaves int // number of leaf slots (== max(inserted symbols, 2) after Complete/LoadTree)
	nodes  []node
	codes  []Code
	minLen uint32

	rootID int32
	cur    int32 // decode cursor; < 0 means "start over from root"
}

// NewTree creates an empty tree; call Restart before Insert.
func NewTree() *Tree {
	return &Tree{}
}

// Restart discards any previous frequencies and sets the capacity (max
// distinct symbols) for the next Insert/Complete cycle.
func (t *Tree) Restart(capacity int) {
	t.capacity = capacity
	t.freqs = t.freqs[:0]
}

// Insert records the next symbol's frequency (symbols are implicitly
// numbered 0, 1, 2... in insertion order). Returns false if capacity is
// already exhausted.
func (t *Tree) Insert(frequency uint32) bool {
	if len(t.freqs) >= t.capacity {
		return false
	}
	t.freqs = append(t.freqs, freqItem{symbol: uint32(len(t.freqs)), freq: frequency})
	return true
}

func bitsForCount(n int) uint {
	b := uint(0)
	for (1 << b) < n {
		b++
	}
	return b
}

// Complete builds the Huffman tree from the frequencies given to Insert and
// returns each inserted symbol's Code, indexed by symbol id. Returns nil if
// nothing was inserted.
func (t *Tree) Complete() []Code {
	inserted := len(t.freqs)
	if inserted == 0 {
		t.codes = nil
		return nil
	}

	n := inserted
	items := append([]freqItem(nil), t.freqs...)
	if n < 2 {
		// See package doc: fresh in-bounds placeholder id, not upstream's
		// out-of-bounds heap slot.
		items = append(items, freqItem{symbol: uint32(n), freq: 0})
		n = 2
	}

	h := &freqHeap{}
	*h = append(*h, items...)
	heap.Init(h)

	t.nodes = make([]node, 2*n-1)
	t.codes = make([]Code, 2*n-1)
	for i := 0; i < n; i++ {
		t.nodes[i] = node{-1, -1}
	}

	// Drop zero-frequency leaves from the combine heap (they still get a
	// node slot and default {0,0} code — a symbol legitimately absent from
	// this particular tree, e.g. a quality value never seen at this read
	// position, is simply never looked up). Mirrors the special-case
	// handling in HuffmanEncoder::Complete.
	if h.Len() == 2 && (*h)[0].freq == 0 {
		(*h)[0].freq = 1
		if (*h)[1].freq == 0 {
			(*h)[1].freq = 1
		}
		heap.Init(h)
	} else {
		for h.Len() > 2 && (*h)[0].freq == 0 {
			heap.Pop(h)
		}
	}

	present := h.Len() // always >= 2: n >= 2 going in, and the loop above never drops below 2

	for i := 0; i < present-1; i++ {
		left := heap.Pop(h).(freqItem)
		right := heap.Pop(h).(freqItem)

		newID := uint32(n + i)
		heap.Push(h, freqItem{symbol: newID, freq: left.freq + right.freq})

		t.nodes[n+i] = node{int32(left.symbol), int32(right.symbol)}
	}

	for i := n + present - 2; i >= n; i-- {
		l, r := t.nodes[i].left, t.nodes[i].right
		t.codes[l] = Code{Code: t.codes[i].Code << 1, Len: t.codes[i].Len + 1}
		t.codes[r] = Code{Code: (t.codes[i].Code << 1) | 1, Len: t.codes[i].Len + 1}
	}

	t.rootID = int32(n + present - 2)
	t.cur = t.rootID
	t.leaves = n

	t.minLen = uint32(n)
	for i := 0; i < n; i++ {
		if t.codes[i].Len > 0 && t.codes[i].Len < t.minLen {
			t.minLen = t.codes[i].Len
		}
	}

	return t.codes[:inserted]
}

// Codes returns the codes computed by the most recent Complete call.
func (t *Tree) Codes() []Code { return t.codes[:min(len(t.codes), t.numRealSymbols())] }

func (t *Tree) numRealSymbols() int { return len(t.freqs) }

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// MinLen returns the shortest code length among inserted symbols, matching
// GetMinLen. (Upstream uses this to size the DecodeFast speedup table,
// which this port omits — see package doc.)
func (t *Tree) MinLen() uint32 { return t.minLen }

// RootID, Leaves, IsLeaf, LeftChild, and RightChild expose the tree shape
// built by Complete for callers that need to walk it directly — e.g.
// realdsrc, which reuses this package's (already bit-exact-verified)
// Complete() for tree construction but serializes the result in real
// DSRC's own on-disk format rather than this package's StoreTree/LoadTree
// shape (see package doc for why those differ).
func (t *Tree) RootID() int32             { return t.rootID }
func (t *Tree) Leaves() int               { return t.leaves }
func (t *Tree) IsLeaf(id int32) bool      { return id < int32(t.leaves) }
func (t *Tree) LeftChild(id int32) int32  { return t.nodes[id].left }
func (t *Tree) RightChild(id int32) int32 { return t.nodes[id].right }

// StoreTree serializes the tree shape as a preorder walk: a 0 bit marks an
// internal node (followed by its left then right subtree), a 1 bit marks a
// leaf (followed by a fixed-width leaf id). Mirrors
// HuffmanEncoder::StoreTree, using this package's single tree shape (see
// package doc) rather than upstream's memory layout, and omitting the
// stored byte-length-prefix (a decode fast-path upstream needs to jump over
// a tree without parsing it, not needed here).
func (t *Tree) StoreTree(w *bitio.Writer) {
	bitsPerID := bitsForCount(t.leaves)
	w.PutWord(uint32(t.leaves))
	w.PutByte(byte(t.minLen))
	t.storeNode(w, t.rootID, bitsPerID)
}

func (t *Tree) storeNode(w *bitio.Writer, id int32, bitsPerID uint) {
	if id < int32(t.leaves) {
		w.PutBit(1)
		w.PutBits(uint32(id), bitsPerID)
		return
	}
	w.PutBit(0)
	t.storeNode(w, t.nodes[id].left, bitsPerID)
	t.storeNode(w, t.nodes[id].right, bitsPerID)
}

// LoadTree reconstructs a Tree from bits written by StoreTree, ready for
// Decode/DecodeSymbol.
func LoadTree(r *bitio.Reader) *Tree {
	t := &Tree{}
	leaves := int(r.GetWord())
	minLen := r.GetByte()
	bitsPerID := bitsForCount(leaves)

	t.leaves = leaves
	t.minLen = uint32(minLen)
	t.nodes = make([]node, 2*leaves-1)
	for i := 0; i < leaves; i++ {
		t.nodes[i] = node{-1, -1}
	}

	next := int32(leaves)
	t.rootID = t.loadNode(r, bitsPerID, &next)
	t.cur = t.rootID
	return t
}

// loadNode reads one subtree in preorder and returns the id it was
// assigned (a leaf's id is the symbol id read directly from the stream; an
// internal node's id is allocated sequentially as it's encountered).
func (t *Tree) loadNode(r *bitio.Reader, bitsPerID uint, next *int32) int32 {
	if r.GetBit() != 0 {
		return int32(r.GetBits(bitsPerID))
	}
	id := *next
	*next++
	left := t.loadNode(r, bitsPerID, next)
	right := t.loadNode(r, bitsPerID, next)
	t.nodes[id] = node{left, right}
	return id
}

// Decode consumes one bit of a code and returns (symbol, true) once
// traversal reaches a leaf, or (-1, false) if more bits are needed. The
// cursor auto-resets to the root at the start of the symbol following one
// that was just decoded, so repeated calls naturally decode a stream of
// symbols. Mirrors HuffmanEncoder::Decode (the single-bit path; see package
// doc re: DecodeFast).
func (t *Tree) Decode(bit uint32) (symbol int32, done bool) {
	if t.cur < 0 {
		t.cur = t.rootID
	}
	if bit != 0 {
		t.cur = t.nodes[t.cur].right
	} else {
		t.cur = t.nodes[t.cur].left
	}
	if t.cur < int32(t.leaves) {
		sym := t.cur
		t.cur = -1
		return sym, true
	}
	return -1, false
}

// DecodeSymbol pulls bits from r via Decode until one full symbol resolves.
func (t *Tree) DecodeSymbol(r *bitio.Reader) int32 {
	for {
		if sym, ok := t.Decode(r.GetBit()); ok {
			return sym
		}
	}
}
