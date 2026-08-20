package realdsrc

import (
	"github.com/fredericlemoine/dsrcgo/internal/bitio"
	"github.com/fredericlemoine/dsrcgo/internal/huffman"
)

// This file ports HuffmanEncoder::StoreTree/LoadTree/Decode's exact on-disk
// format (src/huffman.cpp) — real DSRC's tree serialization, distinct from
// go/internal/huffman's StoreTree/LoadTree, which deliberately uses one
// unified tree shape instead of upstream's two different ones (see that
// package's doc comment). Real DSRC's format needs both:
//
//   - On the encode side, real DSRC's tree-construction algorithm
//     (HuffmanEncoder::Complete) turned out, after empirical verification
//     against a real archive (see realdsrc_test.go), to build
//     structurally identical trees to go/internal/huffman.Tree.Complete()
//     given the same frequencies — same min-heap-by-(frequency,symbol-id)
//     combine order. So encoding reuses that (already tested) Complete(),
//     just walked and serialized in upstream's own wire format via the
//     RootID/Leaves/IsLeaf/LeftChild/RightChild accessors added to that
//     package for this purpose.
//   - On the decode side, real DSRC uses a different, more compact tree
//     representation than go/internal/huffman's: leaf identity is embedded
//     directly as a negative child-pointer value rather than a separate
//     node slot, and internal node ids are handed out via a shared
//     decrementing counter during a preorder walk (HuffmanEncoder::
//     DecodeProcess) rather than the bottom-up combine-order numbering the
//     encode side's tree uses. realHuffmanTree below is a direct,
//     from-scratch translation of that decode-side algorithm — it doesn't
//     reuse go/internal/huffman.Tree at all, since the two representations
//     aren't shaped alike.

// storeRealHuffmanTree serializes t's shape as upstream's HuffmanEncoder::
// StoreTree does: a stored total-length prefix, then root id, symbol
// count, and minimum code length, then a preorder walk of the tree (0 bit
// = internal node followed by its children, 1 bit = leaf followed by a
// fixed-width leaf id), padded to a byte boundary.
//
// Upstream computes the length prefix by writing a placeholder, walking
// the tree, and then seeking back to patch in the real value once known.
// bitio.Writer is append-only, so this walks the tree into a scratch
// buffer first to measure its size, then writes the real header followed
// by that buffer — legwork upstream doesn't need to do because it can
// seek, but the two ways of arriving at the on-disk bytes agree.
//
// The walk can't be measured analytically from n_symbols alone: only
// leaves with positive frequency are ever actually combined into the
// tree (see huffman.Tree.Complete's zero-frequency pruning), so a
// position with a smaller local alphabet than the block-wide symbol
// count produces a smaller tree, even though n_symbols itself — and thus
// bitsPerID, since upstream sizes leaf ids off the full tree_size — stays
// fixed at the block-wide count.
func storeRealHuffmanTree(w *bitio.Writer, t *huffman.Tree) {
	w.AlignByte()

	n := t.Leaves()
	bitsPerID := bitsForCount(n)

	scratch := bitio.NewWriter()
	storeRealNode(scratch, t, t.RootID(), bitsPerID)
	treeBytes := scratch.Bytes()

	memSize := uint32(4+4+4+1) + uint32(len(treeBytes)) // memSize field itself + root_id + n_symbols + min_len + tree bytes

	w.PutWord(memSize)
	w.PutWord(uint32(t.RootID()))
	w.PutWord(uint32(n))
	w.PutByte(byte(t.MinLen()))
	for _, b := range treeBytes {
		w.PutByte(b)
	}
}

func storeRealNode(w *bitio.Writer, t *huffman.Tree, id int32, bitsPerID uint) {
	if t.IsLeaf(id) {
		w.PutBit(1)
		w.PutBits(uint32(id), bitsPerID)
		return
	}
	w.PutBit(0)
	storeRealNode(w, t, t.LeftChild(id), bitsPerID)
	storeRealNode(w, t, t.RightChild(id), bitsPerID)
}

// realTreeNode is upstream's decode-side node shape: left/right are either
// a non-negative index into the same nodes slice (an internal node), or a
// value <= 0 whose negation is a leaf's symbol id.
type realTreeNode struct {
	left, right int32
}

// realHuffmanTree is a decode-only tree built by loadRealHuffmanTree.
type realHuffmanTree struct {
	nodes  []realTreeNode
	root   int32
	minLen uint32
	cur    int32
}

// loadRealHuffmanTree parses bits written by storeRealHuffmanTree (or by
// real dsrc itself), mirroring HuffmanEncoder::LoadTree + DecodeProcess.
func loadRealHuffmanTree(r *bitio.Reader) *realHuffmanTree {
	r.AlignByte()

	_ = r.GetWord() // memSize: upstream uses this to skip the tree without parsing it; not needed here
	origRootID := int32(r.GetWord())
	nSymbols := int(r.GetWord())
	minLen := uint32(r.GetByte())

	bitsPerID := bitsForCount(nSymbols)

	// Upstream remaps the root id for the decode-side node numbering
	// scheme: internal node ids run 1..remappedRoot, handed out by a
	// shared decrementing counter during the preorder walk below.
	remappedRoot := origRootID - int32(nSymbols) + 1

	t := &realHuffmanTree{
		nodes:  make([]realTreeNode, remappedRoot+1),
		root:   remappedRoot,
		minLen: minLen,
		cur:    remappedRoot,
	}

	next := remappedRoot
	decodeRealNode(r, t.nodes, remappedRoot, &next, bitsPerID)
	r.AlignByte()

	return t
}

// decodeRealNode mirrors HuffmanEncoder::DecodeProcess exactly, including
// its shared-counter internal-node id allocation: an id handed to a call
// that turns out to be a leaf is simply never written to nodes[] (leaves
// don't occupy a slot), so the next call — whichever child needs a slot
// next — naturally reuses it.
func decodeRealNode(r *bitio.Reader, nodes []realTreeNode, nodeID int32, next *int32, bitsPerID uint) int32 {
	if r.GetBit() == 0 {
		*next--
		nodes[nodeID].left = decodeRealNode(r, nodes, *next, next, bitsPerID)
		nodes[nodeID].right = decodeRealNode(r, nodes, *next, next, bitsPerID)
		return nodeID
	}
	return -int32(r.GetBits(bitsPerID))
}

// decode consumes one bit and returns (symbol, true) once traversal
// reaches a leaf, mirroring HuffmanEncoder::Decode (the single-bit path;
// upstream's DecodeFast speedup table is a pure performance optimization,
// not implemented here — same simplification go/internal/huffman makes).
func (t *realHuffmanTree) decode(bit uint32) (int32, bool) {
	if t.cur <= 0 {
		t.cur = t.root
	}
	if bit != 0 {
		t.cur = t.nodes[t.cur].right
	} else {
		t.cur = t.nodes[t.cur].left
	}
	if t.cur <= 0 {
		return -t.cur, true
	}
	return -1, false
}

// DecodeSymbol pulls bits from r until one full symbol resolves.
func (t *realHuffmanTree) DecodeSymbol(r *bitio.Reader) int32 {
	for {
		if sym, ok := t.decode(r.GetBit()); ok {
			return sym
		}
	}
}
