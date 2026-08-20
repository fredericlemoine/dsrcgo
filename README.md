# AI-Generated Code Notice
This project was generated using Claude Code as an experiment in automated C++ to Go conversion. While the functional tests pass successfully, the code has not undergone thorough human review and should be used with appropriate scrutiny.

# dsrcgo

A Go port of [DSRC](https://github.com/refresh-bio/DSRC) (DNA Sequence Reads Compressor), built to be **byte-exact interoperable** with the real C++ `dsrc` binary: files this port writes can be decompressed by real `dsrc`, and files real `dsrc` writes can be decompressed by this port. That's the only thing this repository does — there's no separate "own format," no algorithmically-similar-but-incompatible fallback path. Everything here either matches upstream's on-disk bit layout or is scoped out and clearly documented as such.

Coverage is upstream's default compression mode, `-m0` (`-d0 -q0`, no lossy quality, no CRC, no color-space, no field filtering) — see [Scope and missing parts](#scope-and-missing-parts) for exactly what that excludes.

This module is self-contained and has no dependency on the original C++ project at build or run time — the [DSRC C++ source](https://github.com/refresh-bio/DSRC) is only relevant if you want to reproduce the interop verification below yourself (build the real `dsrc` binary and diff its output against this port's), not to build or use `dsrcgo` itself.

## How this was verified

Not by reading the C++ source alone. A checkout of the original DSRC C++ source (`src/BlockCompressor.cpp` and friends) was compiled into a debug build with assertions enabled (no `-DNDEBUG`), which embeds 4-byte position markers between every section of a compressed block (`CONTROL_CHECK_W`/`R`, normally compiled out in release builds). That let section boundaries be located precisely in real archives without first having to implement decoding for later sections — every format decision below was pinned down against real output, not inferred from source reading.

Every package here has tests that decode a real archive (checked into `internal/realdsrc/testdata/`, produced by an actual `dsrc` build) and assert this port's own encoder output matches it exactly, byte for byte, plus decodes it back to confirm round-tripping. Beyond the unit tests, the full CLI has been run end-to-end against a real `dsrc` binary on files up to ~30,000 records / 7 blocks / 7 MB, including records with N/ambiguity codes and variable-length reads — output byte-identical in both directions.

Run `go test ./...` for the full suite (60 tests as of this writing).

## Quick start

```bash
go build -o dsrcgo ./cmd/dsrcgo
go test ./...
```

```bash
./dsrcgo generate demo.fastq 50000 100   # synthetic test data
./dsrcgo pack demo.dsrc demo.fastq       # -> a real dsrc archive
```

To confirm interoperability against a real `dsrc` binary, build one from the [original C++ source](https://github.com/refresh-bio/DSRC) (e.g. `make -f Makefile.osx bin`, producing `bin/dsrc`), then:

```bash
dsrc d demo.dsrc demo.decoded.fastq   # decompressed by the ACTUAL C++ tool
diff demo.fastq demo.decoded.fastq    # identical

dsrc c -m0 demo.fastq demo2.dsrc      # compressed by the ACTUAL C++ tool
./dsrcgo unpack demo2.dsrc demo2.out.fastq
diff demo.fastq demo2.out.fastq       # identical
```

`pack`/`unpack` assume quality offset 33 (Sanger/Illumina 1.8+, essentially all modern data).

## Multi-threading

Blocks are fully independent — encoding or decoding one never depends on another's state (each resets its own DNA/quality/tag model per block, matching upstream's own per-block design). `pack` and `unpack` exploit that: file I/O stays strictly sequential (a single `os.File`, and `archive.Writer.WriteBlock` must be called in block order), but the actual parsing/encoding (`pack`) or decoding (`unpack`) work for every block runs across a worker pool sized to `runtime.NumCPU()` by default, overridable via a trailing `[workers]` argument on either command.

This is a from-scratch application of that independence, not a port of DSRC's own C++ worker/queue pipeline (`src/DsrcWorker.cpp`, `src/DsrcOperator.cpp`), which isn't ported here — a simpler, purely block-parallel scheme was enough to get real, measurable speedup (confirmed via wall-clock vs. CPU-time on a multi-block file — e.g. ~45ms wall time against ~260ms of total CPU time across workers on a 7-block file on a 10-core machine).

`pack` also verifies every block by decoding it back and comparing against the original bytes immediately after encoding, in the same worker goroutine — a subtle interop bug surfaces as an error on that specific block rather than silently producing a corrupt archive.

## Package guide

| Package | Ports | Notes |
|---|---|---|
| [`fastq`](internal/fastq) | `FastqParser`, `IFastqStreamReader` (`src/FastqParser.*`, `src/FastqStream.*`) | Record-boundary-aligned chunk reading and parsing. |
| [`bitio`](internal/bitio) | `BitMemoryWriter`/`Reader` (`src/BitMemory.h`) | MSB-first bit packing. A from-scratch equivalent, not byte-for-byte identical internals — only bit *order* needs to match between a writer and its reader, which this preserves exactly. |
| [`huffman`](internal/huffman) | `HuffmanEncoder`'s tree-construction algorithm (`src/huffman.cpp`'s `Insert`/`Complete`) | Reused by `realdsrc` for its verified min-heap-by-(frequency, symbol-id) tree construction; `realdsrc/realhuffman.go` implements upstream's actual on-disk tree serialization separately, since that differs from any simpler shape this package might otherwise use. |
| [`dna`](internal/dna) | The fixed IUPAC symbol table from `RecordsProcessor.cpp` (A=0, G=1, C=2, T=3, N=4, ...) | Just the table — the base↔index mapping is identical regardless of which DNA scheme (B2 or Huffman) wraps it. |
| [`archive`](internal/archive) | `DsrcFileHeader`/`DsrcFileFooter` (`src/DsrcFile.*`) | The real container format: same magic bytes (`0xAA`/`0xCC`), same version numbers, same field layout — including one genuine quirk found by diffing real archives: the footer's block-size index is raw **native-endian** `uint32`s (a `memcpy`-style write in the C++), while every other multi-byte field in the format is big-endian. |
| [`realdsrc`](internal/realdsrc) | `BlockCompressor` and everything under it (`src/BlockCompressor.cpp`, `src/DnaModeler*`, `src/QualityModelerProxy.h`+`QualityPositionModeler.cpp`+`QualityRLEModeler.cpp`, `src/TagModeler.cpp`) | The actual compression pipeline — see below. |

### `realdsrc` in detail

| File | Covers |
|---|---|
| `realdsrc.go` | `ChunkHeader` (excluding color-space and CRC32 fields). |
| `preprocess.go` | The DNA/quality coupling real DSRC's lossless mode has: an ambiguous base (non-ACGT) at low quality gets dropped from the DNA stream and "smuggled" into the quality stream instead, as a value ≥128 encoding both the base's identity and its quality. |
| `dna.go` | B2 (2-bit ACGT packing) and Huffman DNA schemes, with upstream's exact scheme-selection logic — including one place this port deliberately *doesn't* replicate upstream: real DSRC's B2 selection only checks symbol *count* (≤4), not whether those symbols are actually plain ACGT, so a block with e.g. only A/T/N/one-more (4 symbols, one non-ACGT) makes real DSRC's own 2-bit packing silently truncate and corrupt that base. This port detects that exact case and errors instead. |
| `quality.go`, `quality_truncated.go`, `quality_rle.go` | All three quality schemes (Plain, Truncated, RLE — including RLE's degenerate single-symbol path) and upstream's `rawLength`/`thLength`/`rleLength`-ratio scheme selection. |
| `realhuffman.go` | Upstream's actual `HuffmanEncoder::StoreTree`/decode-tree wire format — a genuinely different on-disk shape from a simpler encode/decode-symmetric tree, reusing `huffman.Tree.Complete()` for construction (see the package guide above) but serializing it upstream's way. |
| `tag.go`, `tag_analyze.go`, `tag_numeric.go`, `tag_text.go` | Tag (FASTQ header/ID) tokenization and per-field modeling: constant fields, numeric fields (`DeltaConst`, `DeltaVar`/`ValueVar` with or without upstream's Huffman-optimized `var_stat_encode` path), and free-text fields (hamming-mask + per-position Huffman). |
| `block.go` | Ties it all together — `EncodeBlock`/`DecodeBlock`, mirroring `BlockCompressor::Store`/`Read`'s exact section order (metadata, tags, quality, DNA) and the quality-length bits `StoreTags` interleaves into the tag section itself. |

## Scope and missing parts

Not implemented, all returning a clear error rather than producing silently-wrong output:

- **`-d>0` / `-q>0` (order-N DNA/quality)** — only the default order-0 schemes are covered.
- **Lossy quality mode (`-l`, Illumina binning)** — `LossyRecordsProcessor` isn't ported.
- **CRC32 (`-c`)** — checksum calculation/verification isn't implemented.
- **Color-space datasets (SOLiD)** — `FLAG_DELTA_CONSTANT` handling isn't ported.
- **Field filtering (`-f`)** — `FastqParserExt`'s tag-field-preservation flags aren't wired through.
- **Tag RLE schemes** (`ValueRle`/`DeltaRle`) — scheme *selection* is fully implemented (needed to correctly detect when real dsrc would pick one and refuse cleanly), only the encoding itself is missing.
- **Mixed-format tag fallback** (`TagRawEncoder`/`Decoder`) — used when a block's tags don't all tokenize to the same field layout.

Also not reproduced, deliberately, because reproducing it would mean matching a bug rather than a format:

- **Huffman's `DecodeFast` speedup table** — a pure performance optimization upstream uses to skip repeated single-bit tree walks; semantically equivalent to (slower) single-bit decoding, which is what's implemented.
- **The single-symbol Huffman edge case's undefined behavior** — upstream's `n_symbols < 2` padding path reads past the end of a size-1 heap allocation in that specific case. This port uses a fresh in-bounds placeholder instead, which round-trips correctly without relying on undefined behavior.

## Repository layout

```
.
├── cmd/dsrcgo/       CLI: pack, unpack, generate
└── internal/
    ├── fastq/         FASTQ chunk reading + parsing
    ├── bitio/          MSB-first bit packing
    ├── huffman/        Huffman tree construction (reused by realdsrc)
    ├── dna/            Fixed IUPAC symbol table (reused by realdsrc)
    ├── archive/         Real DSRC file container format
    └── realdsrc/       The compression pipeline: metadata, DNA, quality, tags, block assembly
```

## License and attribution

This is a derivative work of [DSRC](https://github.com/refresh-bio/DSRC) (© Lucas Roguski and Sebastian Deorowicz), licensed GNU GPL v2 — see [NOTICE](NOTICE) for full attribution and the papers to cite.
