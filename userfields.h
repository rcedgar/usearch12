#ifndef f
#error f not defined
#endif

f(query)			// Query sequence label.
f(target)			// Target sequenc label.
f(clusternr)		// Cluster number.
f(evalue)			// E-value.
f(id)				// Pct identity.
f(fractid)			// Fractional identity.
f(dist)				// Distance (0.0 .. 1.0)
f(mid)				// Fraction of letter pairs that match.
f(pctpv)			// Percent positives (substitution score > 0).
f(pctgaps)			// Percent gap columns.
f(pairs)			// Number of letter pair columns.
f(gaps)				// Number of gapped columns (internal gaps only).
f(allgaps)			// Number of gapped columns (internal + terminal gaps).
f(qlo)				// 0-based start position of alignment in query sequence.
f(qhi)				// 0-based end position of alignment in query sequence.
f(tlo)				// 0-based start position of alignment in target sequence.
f(thi)				// 0-based end position of alignment in target sequence.
f(qlot)				// 0-based start position of alignment in query sequence.
f(qhit)				// 0-based end position of alignment in query sequence.
f(qunt)				// Number of letters unaligned at end of query sequence
f(tlot)				// 0-based start position of alignment in target sequence.
f(thit)				// 0-based end position of alignment in target sequence.
f(tunt)				// Number of letters unaligned at end of query sequence
f(pv)				// Number of positive columns (substitution score > 0).
f(ql)				// Query sequence length.
f(tl)				// Target sequence length.
f(qs)				// Query segment length.
f(ts)				// Target segment length.
f(alnlen)			// Number of alignment columns.
f(opens)			// Number of gap opens.
f(exts)				// Number of gap extensions (does not include gap opens).
f(raw)				// Raw alignment score.
f(bits)				// Bit score.
f(aln)				// Alignment as string M=letter pair, D=delete (gap in query), I=insert (gap in target).
f(caln)				// Compressed alignment (see uc file format).
f(qseq)				// Full-length query sequence.
f(tseq)				// Full-length target sequence.
f(qseg)				// Query segment.
f(tseg)				// Target segment
f(qsegf)			// Query segment with flanks
f(tsegf)			// Target segment with flanks
f(qstrand)			// Query strand.
f(tstrand)			// Target strand.
f(qrow)				// Aligned query segment with gaps.
f(trow)				// Aligned target segment with gaps.
f(qrowdots)			// Aligned query segment with gaps and dots for matched positions.
f(trowdots)			// Aligned target segment with gaps and dots for matched positions.
f(qframe)			// Query frame (-3 to +3).
f(tframe)			// Target frame (-3 to +3).
f(mism)				// Number of mismatches.
f(ids)				// Number of columns with matching letters.
f(qcov)				// Fraction of query sequence covered by alignment.
f(tcov)				// Fraction of target sequence covered by alignment.
f(diffs)			// Number of differences (gaps + mismatches).
f(diffsa)			// Number of differences (gaps + mismatches), alphabetical.
f(editdiffs)		// Number of differences (gaps + mismatches + termgaps).
f(abskew)			// Abundance skew.
f(qlor)				// 0-based query start coordinate.
f(qhir)				// 0-based query end coordinate.
f(tlor)				// 0-based target start coordinate.
f(thir)				// 0-based target end coordinate.
f(orflo)			// Start coordinate of query ORF in nt sequence.
f(orfhi)			// End coordinate of query ORF in nt sequence.
f(orfframe)			// Frame of query ORF.
f(orfseqaa)			// 
f(orfseqnt)			// 
f(orfsegnt)			// 
f(gc)				// 
f(kmerid)				// 
f(qtrimlo)				// 
f(qtrimhi)				// 
f(qtrimseq)				// 

#undef f
