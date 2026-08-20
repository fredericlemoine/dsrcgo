// Command dsrcgo is a real-dsrc-compatible FASTQ compressor: pack writes an
// archive (internal/archive's exact container format wrapping
// internal/realdsrc's exact block format) that the actual C++ dsrc binary
// can decompress; unpack reads an archive from either tool. See ../../README.md
// for exactly which cases are covered and which aren't.
//
// Blocks are fully independent — each resets its own DNA/quality/tag model
// state (see realdsrc.EncodeBlock/DecodeBlock) — so pack and unpack process
// them across a worker pool of goroutines: parsing/encoding (pack) or
// decoding (unpack) happens in parallel, while the actual file write stays
// strictly sequential (archive.Writer.WriteBlock must be called in order,
// and a single os.File can't be written from multiple goroutines safely
// anyway). This is a from-scratch application of that independence, not a
// port of DSRC's own C++ worker/queue pipeline (src/DsrcWorker.cpp,
// src/DsrcOperator.cpp), which isn't ported here.
//
// Usage:
//
//	dsrcgo pack <file.dsrc> <file.fastq> [workers]
//	dsrcgo unpack <file.dsrc> <file.fastq> [workers]
//	dsrcgo generate <file.fastq> <numRecords> <readLen>
package main

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"runtime"
	"strconv"
)

// qualityOffset is assumed throughout: 33 (Sanger/Illumina 1.8+), the value
// real dsrc auto-detects for essentially all modern FASTQ data. See the
// README for what isn't covered (color-space, non-standard offsets, etc).
const qualityOffset = 33

func main() {
	if len(os.Args) < 2 {
		usage()
	}

	var err error
	switch os.Args[1] {
	case "pack":
		if len(os.Args) != 4 && len(os.Args) != 5 {
			usage()
		}
		workers := runtime.NumCPU()
		if len(os.Args) == 5 {
			workers, err = strconv.Atoi(os.Args[4])
		}
		if err == nil {
			err = runPack(os.Args[2], os.Args[3], workers)
		}
	case "unpack":
		if len(os.Args) != 4 && len(os.Args) != 5 {
			usage()
		}
		workers := runtime.NumCPU()
		if len(os.Args) == 5 {
			workers, err = strconv.Atoi(os.Args[4])
		}
		if err == nil {
			err = runUnpack(os.Args[2], os.Args[3], workers)
		}
	case "generate":
		if len(os.Args) != 5 {
			usage()
		}
		var n, readLen int
		n, err = strconv.Atoi(os.Args[3])
		if err == nil {
			readLen, err = strconv.Atoi(os.Args[4])
		}
		if err == nil {
			err = runGenerate(os.Args[2], n, readLen)
		}
	default:
		usage()
	}

	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
}

func usage() {
	fmt.Fprintln(os.Stderr, "usage:")
	fmt.Fprintln(os.Stderr, "  dsrcgo pack <file.dsrc> <file.fastq> [workers]")
	fmt.Fprintln(os.Stderr, "  dsrcgo unpack <file.dsrc> <file.fastq> [workers]")
	fmt.Fprintln(os.Stderr, "  dsrcgo generate <file.fastq> <numRecords> <readLen>")
	os.Exit(2)
}

// runGenerate writes a synthetic but realistic-shaped FASTQ file: an
// Illumina-like incrementing tag, DNA with local k-mer structure, and
// quality scores that decay toward the read's end.
func runGenerate(path string, n, readLen int) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	rng := rand.New(rand.NewSource(42))
	bases := []byte("AGCT")

	for i := 0; i < n; i++ {
		fmt.Fprintf(w, "@SRR001471.%d HWI-EAS7:1:3:%d:%d/1\n", 1000000+i, 1000+rng.Intn(9000), 1000+rng.Intn(9000))

		seq := make([]byte, readLen)
		ctx := 0
		for j := range seq {
			idx := ctx
			if rng.Float64() >= 0.8 {
				idx = rng.Intn(4)
			}
			seq[j] = bases[idx]
			ctx = idx
		}
		w.Write(seq)
		w.WriteByte('\n')
		w.WriteString("+\n")

		qual := make([]byte, readLen)
		cur := 39
		for j := range qual {
			decay := j * 40 / (readLen*3 + 1)
			target := cur
			if target-decay > 0 {
				target -= decay
			}
			cur = target + rng.Intn(3) - 1
			if cur < 0 {
				cur = 0
			}
			if cur > 40 {
				cur = 40
			}
			qual[j] = byte(33 + cur)
		}
		w.Write(qual)
		w.WriteByte('\n')
	}
	if err := w.Flush(); err != nil {
		return err
	}

	fmt.Printf("wrote %d records (read length %d) to %s\n", n, readLen, path)
	return nil
}
