// Package bitio provides MSB-first variable-width bit packing, the
// primitive the Huffman coder (src/huffman.h) needs for its variable-length
// codes and for serializing tree shape. It's a simplified, from-scratch
// equivalent of core::BitMemoryWriter/Reader (src/BitMemory.h): same bit
// order and the same PutBits/GetBits contract, but a single internal
// accumulator instead of upstream's separate 8-bit (reader) and 32-bit
// (writer) word buffers — an implementation detail that doesn't affect the
// resulting bitstream, since only bit order matters for compatibility
// between a Writer and a Reader.
package bitio

// Writer accumulates bits MSB-first and packs them into bytes.
type Writer struct {
	buf   []byte
	acc   uint64
	nbits uint
}

func NewWriter() *Writer {
	return &Writer{}
}

// PutBits writes the low n bits of word, most significant first. n must be
// in [0, 32]; n == 0 is a no-op.
func (w *Writer) PutBits(word uint32, n uint) {
	if n == 0 {
		return
	}
	word &= uint32((uint64(1) << n) - 1)
	w.acc = (w.acc << n) | uint64(word)
	w.nbits += n
	for w.nbits >= 8 {
		w.nbits -= 8
		w.buf = append(w.buf, byte(w.acc>>w.nbits))
	}
}

func (w *Writer) PutBit(bit uint32) { w.PutBits(bit&1, 1) }
func (w *Writer) PutByte(b byte)    { w.PutBits(uint32(b), 8) }
func (w *Writer) PutWord(v uint32)  { w.PutBits(v, 32) }

// WriteByte satisfies io.ByteWriter (and rangecoder.ByteWriter), so a
// Writer can be handed directly to a range coder — bytes it writes are
// packed at whatever bit position the writer is currently at, same as any
// other PutBits call, so range-coded and bit-packed data can share one
// stream with no explicit byte-alignment step between them.
func (w *Writer) WriteByte(b byte) error {
	w.PutByte(b)
	return nil
}

// AlignByte pads any pending partial byte with zero bits, flushing it
// immediately, mid-stream — mirrors BitMemoryWriter::FlushPartialWordBuffer
// (src/BitMemory.h), which several DSRC formats call between sub-sections
// of a bit-packed stream to force byte alignment. A no-op if already
// byte-aligned.
func (w *Writer) AlignByte() {
	if w.nbits > 0 {
		w.buf = append(w.buf, byte(w.acc<<(8-w.nbits)))
		w.acc = 0
		w.nbits = 0
	}
}

// Bytes flushes any partial trailing byte (zero-padded in the low bits) and
// returns the packed data.
func (w *Writer) Bytes() []byte {
	if w.nbits > 0 {
		w.buf = append(w.buf, byte(w.acc<<(8-w.nbits)))
		w.acc = 0
		w.nbits = 0
	}
	return w.buf
}

// Reader unpacks bits MSB-first from a byte slice, mirroring Writer. Reading
// past the end of data yields zero bits, matching a Writer's zero padding.
type Reader struct {
	data  []byte
	pos   int
	acc   uint64
	nbits uint
}

func NewReader(data []byte) *Reader {
	return &Reader{data: data}
}

func (r *Reader) fill(n uint) {
	for r.nbits < n {
		var b byte
		if r.pos < len(r.data) {
			b = r.data[r.pos]
			r.pos++
		}
		r.acc = (r.acc << 8) | uint64(b)
		r.nbits += 8
	}
}

// GetBits reads n bits (n in [0, 32]) and returns them right-aligned.
func (r *Reader) GetBits(n uint) uint32 {
	if n == 0 {
		return 0
	}
	r.fill(n)
	r.nbits -= n
	return uint32((r.acc >> r.nbits) & ((uint64(1) << n) - 1))
}

func (r *Reader) GetBit() uint32 { return r.GetBits(1) }
func (r *Reader) GetByte() byte  { return byte(r.GetBits(8)) }

// ReadByte satisfies io.ByteReader (and rangecoder.ByteReader) — see
// Writer.WriteByte.
func (r *Reader) ReadByte() (byte, error) {
	return r.GetByte(), nil
}

// AlignByte discards any pending partial-byte bits, moving the read cursor
// to the next byte boundary — mirrors BitMemoryReader::FlushInputWordBuffer.
// A no-op if already byte-aligned.
func (r *Reader) AlignByte() {
	r.nbits = 0
}
func (r *Reader) GetWord() uint32 { return r.GetBits(32) }
