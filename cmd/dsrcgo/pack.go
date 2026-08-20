package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"sync"
	"time"

	"github.com/fredericlemoine/dsrcgo/internal/archive"
	"github.com/fredericlemoine/dsrcgo/internal/fastq"
	"github.com/fredericlemoine/dsrcgo/internal/realdsrc"
)

// runPack writes a real-dsrc-compatible archive from a FASTQ file. Chunks
// are read from disk sequentially (fastq.ChunkReader reuses its internal
// buffer, so it isn't safe for concurrent reads), then parsed, encoded, and
// immediately decoded back for verification across a worker pool — only
// the final archive.Writer.WriteBlock calls are sequential, preserving
// block order.
func runPack(archivePath, fastqPath string, workers int) error {
	if workers < 1 {
		workers = 1
	}

	chunks, rawSize, err := readChunks(fastqPath)
	if err != nil {
		return err
	}
	if len(chunks) == 0 {
		return fmt.Errorf("pack: %s has no records", fastqPath)
	}

	results := make([]encodeResult, len(chunks))

	start := time.Now()

	sem := make(chan struct{}, workers)
	var wg sync.WaitGroup
	for i, c := range chunks {
		wg.Add(1)
		sem <- struct{}{}
		go func(i int, c []byte) {
			defer wg.Done()
			defer func() { <-sem }()
			results[i] = encodeAndVerifyChunk(c)
		}(i, c)
	}
	wg.Wait()

	w, err := archive.Create(archivePath)
	if err != nil {
		return err
	}
	w.SetDatasetType(archive.DatasetType{QualityOffset: qualityOffset})
	w.SetCompressionSettings(archive.CompressionSettings{})

	var records int
	for i, r := range results {
		if r.err != nil {
			w.Close()
			os.Remove(archivePath)
			return fmt.Errorf("block %d: %w", i, r.err)
		}
		if err := w.WriteBlock(r.data, r.numRecords); err != nil {
			w.Close()
			os.Remove(archivePath)
			return fmt.Errorf("writing block %d: %w", i, err)
		}
		records += r.numRecords
	}
	if err := w.Close(); err != nil {
		return err
	}

	elapsed := time.Since(start)
	info, err := os.Stat(archivePath)
	if err != nil {
		return err
	}

	fmt.Printf("%s -> %s\n", fastqPath, archivePath)
	fmt.Printf("  blocks:       %d\n", len(chunks))
	fmt.Printf("  records:      %d\n", records)
	fmt.Printf("  original:     %d bytes\n", rawSize)
	fmt.Printf("  archive:      %d bytes (%.2fx)\n", info.Size(), float64(rawSize)/float64(info.Size()))
	fmt.Printf("  workers:      %d\n", workers)
	fmt.Printf("  wall time:    %s\n", elapsed)
	fmt.Println("  round-trip verified against real dsrc's own bit layout: OK")
	return nil
}

type encodeResult struct {
	data       []byte
	numRecords int
	err        error
}

// encodeAndVerifyChunk encodes one chunk, then immediately decodes it back
// and compares against the original bytes — the same defensive pattern as
// every real-archive test in internal/realdsrc, applied here so a subtle
// interop bug surfaces immediately rather than producing a corrupt archive.
func encodeAndVerifyChunk(chunk []byte) (r encodeResult) {
	records, _, err := fastq.ParseChunk(chunk)
	if err != nil {
		r.err = fmt.Errorf("parsing: %w", err)
		return r
	}

	data, err := realdsrc.EncodeBlock(chunk, qualityOffset)
	if err != nil {
		r.err = fmt.Errorf("encoding: %w", err)
		return r
	}

	decoded, err := realdsrc.DecodeBlock(data, qualityOffset, false)
	if err != nil {
		r.err = fmt.Errorf("decoding for verification: %w", err)
		return r
	}
	// DecodeBlock always terminates the last record with a newline;
	// ChunkReader trims exactly that one trailing byte from its chunks.
	if !bytes.Equal(decoded, append(append([]byte(nil), chunk...), '\n')) {
		r.err = fmt.Errorf("round-trip mismatch")
		return r
	}

	return encodeResult{data: data, numRecords: len(records)}
}

// readChunks splits the file into record-boundary-aligned chunks up front.
func readChunks(path string) (chunks [][]byte, rawSize int64, err error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, 0, err
	}
	defer f.Close()

	cr := fastq.NewChunkReader(f, fastq.DefaultChunkSize)
	for {
		c, err := cr.ReadNextChunk()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, 0, fmt.Errorf("reading chunk %d: %w", len(chunks), err)
		}
		rawSize += int64(len(c))
		chunks = append(chunks, append([]byte(nil), c...))
	}
	return chunks, rawSize, nil
}

// runUnpack reads a dsrc-format archive (from real dsrc, or from runPack)
// and writes the reconstructed FASTQ text. Blocks are read from the
// archive sequentially (single os.File), then decoded across a worker
// pool, then written to the output file in order.
func runUnpack(archivePath, fastqPath string, workers int) error {
	if workers < 1 {
		workers = 1
	}

	r, err := archive.Open(archivePath)
	if err != nil {
		return err
	}

	var blocks [][]byte
	for {
		data, err := r.ReadNextBlock()
		if err == io.EOF {
			break
		}
		if err != nil {
			r.Close()
			return fmt.Errorf("reading block %d: %w", len(blocks), err)
		}
		blocks = append(blocks, data)
	}
	r.Close()

	type result struct {
		chunk []byte
		err   error
	}
	results := make([]result, len(blocks))

	start := time.Now()

	sem := make(chan struct{}, workers)
	var wg sync.WaitGroup
	for i, b := range blocks {
		wg.Add(1)
		sem <- struct{}{}
		go func(i int, b []byte) {
			defer wg.Done()
			defer func() { <-sem }()
			chunk, err := realdsrc.DecodeBlock(b, qualityOffset, false)
			results[i] = result{chunk: chunk, err: err}
		}(i, b)
	}
	wg.Wait()

	out, err := os.Create(fastqPath)
	if err != nil {
		return err
	}
	defer out.Close()

	w := bufio.NewWriter(out)
	for i, res := range results {
		if res.err != nil {
			return fmt.Errorf("decoding block %d: %w", i, res.err)
		}
		if _, err := w.Write(res.chunk); err != nil {
			return err
		}
	}
	if err := w.Flush(); err != nil {
		return err
	}

	fmt.Printf("%s -> %s\n", archivePath, fastqPath)
	fmt.Printf("  blocks:    %d\n", len(blocks))
	fmt.Printf("  workers:   %d\n", workers)
	fmt.Printf("  wall time: %s\n", time.Since(start))
	return nil
}
