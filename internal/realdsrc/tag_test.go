package realdsrc

import (
	"bytes"
	"fmt"
	"math/rand"
	"testing"

	"github.com/fredericlemoine/dsrcgo/internal/bitio"
)

func roundTripTags(t *testing.T, tags [][]byte) []byte {
	t.Helper()
	// These tests exercise tags in isolation with a fixed quality length
	// (min==max), so the quality-length interleaving EncodeTags/DecodeTags
	// do for real blocks (see tag.go's doc comment) never actually reads
	// or writes any bits here.
	w := bitio.NewWriter()
	if err := EncodeTags(w, tags, nil, 0, 0); err != nil {
		t.Fatalf("EncodeTags: %v", err)
	}
	data := w.Bytes()

	got, _, err := DecodeTags(bitio.NewReader(data), len(tags), 0, 0)
	if err != nil {
		t.Fatalf("DecodeTags: %v", err)
	}
	if len(got) != len(tags) {
		t.Fatalf("got %d tags, want %d", len(got), len(tags))
	}
	for i := range tags {
		if !bytes.Equal(got[i], tags[i]) {
			t.Fatalf("tag %d mismatch\n got  %q\n want %q", i, got[i], tags[i])
		}
	}
	return data
}

func TestTagsAllConstant(t *testing.T) {
	tags := make([][]byte, 20)
	for i := range tags {
		tags[i] = []byte("@instrument:1:1:1000:2000")
	}
	roundTripTags(t, tags)
}

func TestTagsDeltaConstCounter(t *testing.T) {
	tags := make([][]byte, 500)
	for i := range tags {
		tags[i] = []byte(fmt.Sprintf("@READ.%d", 1000000+i))
	}
	roundTripTags(t, tags)
}

func TestTagsDeltaVarNonConstantStep(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	tags := make([][]byte, 300)
	val := 1000000
	for i := range tags {
		val += 1 + rng.Intn(9)
		tags[i] = []byte(fmt.Sprintf("@READ.%d", val))
	}
	roundTripTags(t, tags)
}

func TestTagsValueVarWideRange(t *testing.T) {
	rng := rand.New(rand.NewSource(2))
	tags := make([][]byte, 300)
	for i := range tags {
		v := 100000 + rng.Intn(900000) // wide, non-monotonic range
		tags[i] = []byte(fmt.Sprintf("@READ.%d", v))
	}
	roundTripTags(t, tags)
}

func TestTagsValueVarSmallCardinality(t *testing.T) {
	// Few distinct values within a wide numeric range: should trigger
	// real dsrc's var_stat_encode (Huffman-optimized) path.
	rng := rand.New(rand.NewSource(3))
	values := []int{100000, 250000, 500000, 750000, 999999}
	tags := make([][]byte, 300)
	for i := range tags {
		v := values[rng.Intn(len(values))]
		tags[i] = []byte(fmt.Sprintf("@READ.%d", v))
	}
	roundTripTags(t, tags)
}

func TestTagsTextFieldVaryingChars(t *testing.T) {
	rng := rand.New(rand.NewSource(4))
	letters := []byte("ACGTN")
	tags := make([][]byte, 200)
	for i := range tags {
		bc := make([]byte, 8)
		for j := range bc {
			bc[j] = letters[rng.Intn(len(letters))]
		}
		tags[i] = append([]byte("@BC."), bc...)
	}
	roundTripTags(t, tags)
}

func TestTagsIlluminaLike(t *testing.T) {
	rng := rand.New(rand.NewSource(5))
	tags := make([][]byte, 1000)
	for i := range tags {
		x := 1000 + rng.Intn(9000)
		y := 1000 + rng.Intn(9000)
		tags[i] = []byte(fmt.Sprintf("@SRR001471.%d HWI-EAS7:1:3:%d:%d/1", 1000000+i, x, y))
	}
	roundTripTags(t, tags)
}

func TestTagsMixedFormattingErrors(t *testing.T) {
	tags := [][]byte{
		[]byte("@a:b:c"),
		[]byte("@a:b"),
	}
	w := bitio.NewWriter()
	if err := EncodeTags(w, tags, nil, 0, 0); err == nil {
		t.Fatal("expected an error for mixed-formatting tags")
	}
}

func TestTagsSingleTag(t *testing.T) {
	roundTripTags(t, [][]byte{[]byte("@ONLY_ONE lane1:2:3")})
}

func TestTagsVariableLengthTextField(t *testing.T) {
	rng := rand.New(rand.NewSource(6))
	tags := make([][]byte, 200)
	for i := range tags {
		n := 3 + rng.Intn(10)
		name := make([]byte, n)
		for j := range name {
			name[j] = byte('a' + rng.Intn(26))
		}
		tags[i] = append([]byte("@N."), name...)
	}
	roundTripTags(t, tags)
}

func TestTagsLongTextFieldPastCatchAllBoundary(t *testing.T) {
	// A text field longer than maxFieldStatLen (128) exercises the
	// catch-all Huffman bucket.
	rng := rand.New(rand.NewSource(9))
	tags := make([][]byte, 50)
	for i := range tags {
		n := 150
		s := make([]byte, n)
		for j := range s {
			s[j] = byte('a' + rng.Intn(4))
		}
		tags[i] = append([]byte("@L."), s...)
	}
	roundTripTags(t, tags)
}
