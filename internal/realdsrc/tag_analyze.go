// Tag field analysis, ported from TagAnalyzer
// (src/TagModeler.cpp:InitializeFieldsStats/UpdateFieldsStats/
// UpdateNumericField/FinalizeFieldsStats) — the pass BlockCompressor::
// AnalyzeTags runs over every record in a block before StoreTags, deciding
// per field whether it's constant, numeric (and which of 5 numeric
// sub-schemes), or free text.
//
// This is a close, function-by-function translation rather than a
// restructured one, including two upstream quirks that are easy to miss
// but change the actual numbers: BlockCompressor::AnalyzeTags calls
// InitializeFieldsStats(records[0]) and then loops UpdateFieldsStats over
// EVERY record including records[0] again — so the first record's value
// gets double-counted into numValues/deltaValues (only a Huffman-frequency
// nicety, but it changes the exact codes if the Huffman-optimized numeric
// scheme ends up selected) — and UpdateNumericField's delta tracking
// starts at the *second* UpdateFieldsStats call (recordCounter==1, i.e.
// real record[1] vs record[0]), not the first, so the redundant pass over
// record[0] only serves to seed prevValue correctly.
package realdsrc

import (
	"bytes"
)

// separator set, mirroring TagAnalyzer's " ._,=:/-#".
var tagSeparator = func() [256]bool {
	var m [256]bool
	for _, c := range []byte(" ._,=:/-#") {
		m[c] = true
	}
	return m
}()

func splitTagFields(tag []byte) (fields [][]byte, seps []byte) {
	start := 0
	for i := 0; i <= len(tag); i++ {
		end := i == len(tag)
		if !end && !tagSeparator[tag[i]] {
			continue
		}
		fields = append(fields, tag[start:i])
		if !end {
			seps = append(seps, tag[i])
		}
		start = i + 1
	}
	return fields, seps
}

// isNum mirrors core::is_num: a canonical (no leading zero unless "0")
// unsigned decimal integer.
func isNum(s []byte) (int32, bool) {
	if len(s) == 0 {
		return 0, false
	}
	if len(s) > 1 && s[0] == '0' {
		return 0, false
	}
	var v int64
	for _, c := range s {
		if c < '0' || c > '9' {
			return 0, false
		}
		v = v*10 + int64(c-'0')
		if v > 1<<30 {
			return 0, false // outside the range Field's int32 min/max sentinels (±1<<30) can represent
		}
	}
	return int32(v), true
}

const maxFieldStatLen = 128 // Field::MAX_FIELD_STAT_LEN
const maxNumValHuf = 512    // Field::MAX_NUM_VAL_HUF

// Numeric scheme ids, matching Field::NumericSchemeEnum exactly (including
// starting at None=0, which is never actually stored).
const (
	numSchemeNone = iota
	numSchemeValueVar
	numSchemeValueRle
	numSchemeDeltaVar
	numSchemeDeltaRle
	numSchemeDeltaConst
)

type tagField struct {
	sep  byte
	data []byte // first record's raw bytes for this field
	len  int    // len(data)

	isConstant     bool
	isLenConstant  bool
	minLen, maxLen int

	isNumeric          bool
	minValue, maxValue int32
	minDelta, maxDelta int32

	hamMask  []bool
	charFreq [][256]uint32 // index 0..min(maxLen,maxFieldStatLen); position maxFieldStatLen is the catch-all bucket

	// running RLE state during analysis (also replayed during encode)
	rleValCurSym, rleValCurLen, rleValRunLen       int32
	rleDeltaCurSym, rleDeltaCurLen, rleDeltaRunLen int32

	numValues   map[int32]uint32
	deltaValues map[int32]uint32

	// finalized decisions
	isDeltaCoding    bool
	isDeltaConst     bool
	tryRleVal        bool
	tryRleDelta      bool
	numericScheme    byte
	noOfBitsPerNum   uint
	noOfBitsPerValue uint
	varStatEncode    bool
	noOfBitsPerLen   uint
}

type tagAnalysis struct {
	mixedFormatting bool
	fields          []*tagField
	recordCount     int32
}

// analyzeTags runs TagAnalyzer's full pass over a block's tags, mirroring
// BlockCompressor::AnalyzeTags -> TagAnalyzer::{Initialize,Update,Finalize}FieldsStats.
func analyzeTags(tags [][]byte) *tagAnalysis {
	if mixedFormat(tags) {
		return &tagAnalysis{mixedFormatting: true}
	}

	a := initFieldsStats(tags[0])

	prevValues := make([]int32, len(a.fields))
	var recordCounter int32
	for _, tag := range tags { // includes tags[0] again, matching upstream's actual loop
		fields, _ := splitTagFields(tag) // mixedFormat already confirmed every tag tokenizes identically
		for c, f := range fields {
			updateFieldStats(a.fields[c], f, recordCounter, prevValues[c])
			prevValues[c] = fieldValueForUpdate(a.fields[c], f, prevValues[c])
		}
		recordCounter++
	}
	a.recordCount = recordCounter

	finalizeFieldsStats(a)
	return a
}

// mixedFormat reports whether every tag doesn't tokenize identically
// (same field count and separators) — upstream discovers this lazily
// during the analysis loop itself; this package checks it up front to
// keep the loop simpler, with the same outcome.
func mixedFormat(tags [][]byte) bool {
	fields0, seps0 := splitTagFields(tags[0])
	for _, t := range tags[1:] {
		fields, seps := splitTagFields(t)
		if len(fields) != len(fields0) || !bytes.Equal(seps, seps0) {
			return true
		}
	}
	return false
}

// initFieldsStats mirrors TagAnalyzer::InitializeFieldsStats.
func initFieldsStats(tag []byte) *tagAnalysis {
	rawFields, seps := splitTagFields(tag)
	fields := make([]*tagField, len(rawFields))
	for i, raw := range rawFields {
		f := &tagField{
			data:          append([]byte(nil), raw...),
			len:           len(raw),
			minLen:        len(raw),
			maxLen:        len(raw),
			isConstant:    true,
			isLenConstant: true,
			hamMask:       make([]bool, len(raw)),
		}
		for j := range f.hamMask {
			f.hamMask[j] = true
		}
		if i < len(seps) {
			f.sep = seps[i]
		} else {
			// The last field has no real separator; upstream's sep byte
			// for it is whatever raw buffer byte happens to follow the
			// tag text in FastqParser's memory-mapped chunk — the line's
			// own terminator, '\n' (0x0a) for the overwhelming common
			// case of Unix line endings, confirmed empirically against a
			// real archive (a single-field tag "@R1" stored sep=0x0a).
			f.sep = '\n'
		}
		if v, ok := isNum(raw); ok {
			f.isNumeric = true
			f.minValue, f.maxValue = v, v
			f.numValues = map[int32]uint32{v: 1}
			f.minDelta, f.maxDelta = 1<<30, -(1 << 30)
			f.deltaValues = map[int32]uint32{}
		}
		fields[i] = f
	}
	return &tagAnalysis{fields: fields}
}

// fieldValueForUpdate re-derives the numeric value just parsed (if any),
// used as the next call's prevValue — mirrors UpdateFieldsStats setting
// prevFieldValues[c_field] = value only when is_numeric.
func fieldValueForUpdate(f *tagField, raw []byte, prev int32) int32 {
	if f.isNumeric {
		if v, ok := isNum(raw); ok {
			return v
		}
	}
	return prev
}

// updateFieldStats mirrors TagAnalyzer::UpdateFieldsStats's per-field body
// (length/char/constant/ham-mask bookkeeping) plus UpdateNumericField.
func updateFieldStats(f *tagField, raw []byte, recordCounter int32, prevValue int32) {
	n := len(raw)

	if n > f.maxLen {
		f.maxLen = n
	} else if n < f.minLen {
		f.minLen = n
	}

	if len(f.charFreq) < f.maxLen {
		grown := make([][256]uint32, f.maxLen)
		copy(grown, f.charFreq)
		f.charFreq = grown
	}
	charsLen := n
	if charsLen > maxFieldStatLen {
		charsLen = maxFieldStatLen
	}
	for x := 0; x < charsLen; x++ {
		f.charFreq[x][raw[x]]++
	}
	if n > maxFieldStatLen {
		if len(f.charFreq) <= maxFieldStatLen {
			grown := make([][256]uint32, maxFieldStatLen+1)
			copy(grown, f.charFreq)
			f.charFreq = grown
		}
		for x := maxFieldStatLen; x < n; x++ {
			f.charFreq[maxFieldStatLen][raw[x]]++
		}
	}

	if f.isConstant {
		if n != f.len || !bytes.Equal(f.data, raw) {
			f.isConstant = false
		}
	}
	if f.isLenConstant {
		f.isLenConstant = f.len == n
	}

	if f.isNumeric {
		if v, ok := isNum(raw); ok {
			updateNumericField(f, v, prevValue, recordCounter)
		} else {
			f.isNumeric = false
		}
	}

	if !f.isConstant {
		limit := n
		if f.len < limit {
			limit = f.len
		}
		for p := 0; p < limit; p++ {
			f.hamMask[p] = f.hamMask[p] && f.data[p] == raw[p]
		}
	}
}

func updateNumericField(f *tagField, curValue, prevValue, recordCounter int32) {
	if curValue < f.minValue {
		f.minValue = curValue
	} else if curValue > f.maxValue {
		f.maxValue = curValue
	}

	if recordCounter > 0 {
		if f.rleValCurSym != curValue {
			f.rleValRunLen++
			f.rleValCurSym = curValue
			f.rleValCurLen = 0
		} else {
			f.rleValCurLen++
			if f.rleValCurLen > 255 {
				f.rleValCurLen = 0
				f.rleValRunLen++
			}
		}
		if len(f.numValues) > 0 {
			f.numValues[curValue]++
			if len(f.numValues) > maxNumValHuf {
				f.numValues = map[int32]uint32{}
			}
		}
	} else {
		f.rleValCurSym = curValue
		f.rleValCurLen = 0
		f.rleValRunLen = 0
		f.numValues[curValue]++
	}

	if recordCounter >= 1 {
		dvalue := curValue - prevValue
		if recordCounter > 1 {
			if dvalue > f.maxDelta {
				f.maxDelta = dvalue
			} else if dvalue < f.minDelta {
				f.minDelta = dvalue
			}

			if f.rleDeltaCurSym != dvalue {
				f.rleDeltaRunLen++
				f.rleDeltaCurSym = dvalue
				f.rleDeltaCurLen = 0
			} else {
				f.rleDeltaCurLen++
				if f.rleDeltaCurLen > 255 {
					f.rleDeltaCurLen = 0
					f.rleDeltaRunLen++
				}
			}

			if len(f.deltaValues) > 0 {
				f.deltaValues[dvalue]++
				if len(f.deltaValues) > maxNumValHuf {
					f.deltaValues = map[int32]uint32{}
				}
			}
		} else {
			f.maxDelta = dvalue
			f.minDelta = dvalue
			f.rleDeltaCurSym = dvalue
			f.rleDeltaCurLen = 0
			f.rleDeltaRunLen = 0
			f.deltaValues[dvalue]++
		}
	}
}

// finalizeFieldsStats mirrors TagAnalyzer::FinalizeFieldsStats.
func finalizeFieldsStats(a *tagAnalysis) {
	for _, f := range a.fields {
		if !f.isNumeric {
			if !f.isConstant {
				capLen := f.maxLen
				if capLen > maxFieldStatLen {
					capLen = maxFieldStatLen + 1
				}
				if len(f.charFreq) < capLen {
					grown := make([][256]uint32, capLen)
					copy(grown, f.charFreq)
					f.charFreq = grown
				}
				f.noOfBitsPerLen = bitLength(uint32(f.maxLen - f.minLen))
			}
			continue
		}

		var diff int32
		if f.maxValue-f.minValue < f.maxDelta-f.minDelta {
			f.isDeltaCoding = false
			diff = f.maxValue - f.minValue
		} else {
			f.isDeltaCoding = true
			diff = f.maxDelta - f.minDelta
		}

		// Final RLE run-count fixup, mirroring the trailing
		// `lens.push_back(cur_len); if (cur_len>0) { cur_len=0; run_len++; }`:
		// run_len tracks transitions (a run's *start*), so it's missing the
		// very first run only when the stream never transitioned away from
		// it — that first run began via initialization, not a mid-loop
		// transition, so nothing ever incremented run_len for it — UNLESS
		// the trailing (currently-open) run itself has repeats (cur_len>0),
		// in which case run_len already reflects it correctly and no fixup
		// is needed. This does not track the true run count in general
		// (only used here for the ratio heuristic below, matching upstream
		// exactly including this quirk).
		valRuns := f.rleValRunLen
		if f.rleValCurLen > 0 {
			valRuns++
		}
		f.tryRleVal = valRuns > 0 && float64(a.recordCount)/float64(valRuns) > 1.25

		if f.isDeltaCoding {
			f.isDeltaConst = diff == 0
			if !f.isDeltaConst {
				deltaRuns := f.rleDeltaRunLen
				if f.rleDeltaCurLen > 0 {
					deltaRuns++
				}
				f.tryRleDelta = deltaRuns > 0 && float64(a.recordCount)/float64(deltaRuns) > 1.25
			}
		}

		switch {
		case f.isDeltaCoding && f.isDeltaConst:
			f.numericScheme = numSchemeDeltaConst
		case f.isDeltaCoding && f.tryRleDelta:
			f.numericScheme = numSchemeDeltaRle
		case f.tryRleVal:
			f.numericScheme = numSchemeValueRle
		case f.isDeltaCoding:
			f.numericScheme = numSchemeDeltaVar
			d := uint32(f.maxDelta-f.minDelta) + 1
			f.varStatEncode = d <= maxNumValHuf && len(f.deltaValues) > 0
		default:
			f.numericScheme = numSchemeValueVar
			d := uint32(f.maxValue-f.minValue) + 1
			f.varStatEncode = d <= maxNumValHuf && len(f.numValues) > 0
		}

		f.noOfBitsPerNum = bitLength(uint32(diff))
		f.noOfBitsPerValue = bitLength(uint32(f.maxValue - f.minValue))
	}
}
