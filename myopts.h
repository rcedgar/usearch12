#ifndef MY_VERSION
#define MY_VERSION	"12.0"
#endif

#define PROGRAM_NAME	"usearch"

#define A(x)		STR_OPT(x)
#include "cmds.h"

#ifndef VECTOR_OPT
#define VECTOR_OPT(x)	/* empty */
#endif

STR_OPT(log)			// Log file for informational messages.
STR_OPT(userout)		// Hits in tabbed text as specified by -userfields.
STR_OPT(uc)				// Top hit/centroid assignment in uc format.
STR_OPT(alnout)			// Hits as human-readable alignments.
STR_OPT(blast6out)		// Hits in -outfmt 6 BLAST+ output (same as -m8 for older BLAST).
STR_OPT(fastapairs)		// Hits as pair-wise alignments in FASTA format.
STR_OPT(qsegout)		// Query sequence hit segments in FASTA format.
STR_OPT(tsegout)		// Target sequence hit segments in FASTA format.

STR_OPT(output)			// Output filename
STR_OPT(output2)		// Output filename
STR_OPT(report)			// Output filename
STR_OPT(centroids)		// Cluster centroid sequences in FASTA format.
STR_OPT(centroids_fastq)		// Cluster centroid sequences in FASTQ format.
STR_OPT(fastaout)		// FASTA output from some commands.
STR_OPT(matched)		// Query sequences that matched the database.
STR_OPT(matchedfq)		// Query sequences that matched the database.
STR_OPT(notmatched)		// Query sequences that did not match the database.
STR_OPT(notmatchedfq)	// Query sequences that did not match the database.
STR_OPT(dbmatched)		// Database sequences that matched at least one query.
STR_OPT(dbnotmatched)	// Database sequences that did not match a query.
STR_OPT(uchimeout)		// UCHIME output.
STR_OPT(uchimealnout)	// UCHIME alignments.
STR_OPT(uparsealnout)	// UPARSE alignment output.
STR_OPT(otus)			// OTU representative sequences
STR_OPT(zotus)			// OTU representative sequences
STR_OPT(chimeras)		// Sequences classified as chimeric in FASTx format.
STR_OPT(nonchimeras)	// Sequences classified as not chimeric in FASTx format.
//STR_OPT(format)			// Format for sparse distance matrix.
STR_OPT(fastqout_notmerged_fwd)		// Unmerged fwd reads (FASTQ format).
STR_OPT(fastqout_notmerged_rev)		// Unmerged rev reads (FASTQ format).
STR_OPT(fastaout_notmerged_fwd)		// Unmerged fwd reads (FASTA format).
STR_OPT(fastaout_notmerged_rev)		// Unmerged rev reads (FASTA format).
STR_OPT(fastaout_overlap_fwd)		// Aligned seg of fwd reads (FASTA format).
STR_OPT(fastqout_overlap_fwd)		// Aligned seg of fwd reads (FASTQ format).
STR_OPT(fastaout_overlap_rev)		// Aligned seg of fwd reads (FASTA format).
STR_OPT(fastqout_overlap_rev)		// Aligned seg of fwd reads (FASTQ format).
STR_OPT(fastaout_discarded)		// Discarded reads (FASTA format).
STR_OPT(fastqout_discarded)		// Discarded reads (FASTQ format).
STR_OPT(clusters)		// Output filename prefix, one FASTA per cluster.
STR_OPT(tabbedout)		//
STR_OPT(ampout)			//
STR_OPT(alpha)			// String specifying compressed alphabet
STR_OPT(otutabin)		// 
STR_OPT(otutabout)		// 
STR_OPT(biomout)		// 
STR_OPT(hitsout)
STR_OPT(start_motif)
STR_OPT(end_motif)
STR_OPT(fragout)
STR_OPT(dbotus)
STR_OPT(dataotus)
STR_OPT(truncstr)		// truncate FASTA label at first occurence of this string

STR_OPT(matrix)			// Input filename, substitution matrix in BLAST format.
STR_OPT(gapopen)		// String specifying gap open penalty(ies).
STR_OPT(gapext)			// String specifying gap extension penalty(ies).

STR_OPT(strand)			// Strand for nt search (plus or both).
STR_OPT(db)				// Database filename
STR_OPT(input)			// -
STR_OPT(reverse)		// Input filename, reverse reads for merging.
STR_OPT(dbmask)			// Database masking (none, seg, dust, fastamino, fastnucleo, soft).
STR_OPT(userfields)		// Fields to appear in -userout file.
STR_OPT(pattern)		// Pattern for spaced alignment seeds.
STR_OPT(relabel)		// Re-label input sequences in output files.
STR_OPT(fastqout)		// Output file, FASTQ format.
STR_OPT(eetabbedout)	// Output file, expected errors in tabbed text format.
STR_OPT(xdrop_save)		// -
STR_OPT(label_suffix)			// Suffix appended to label of merged reads.
STR_OPT(join_padgap)			// Letters to use for FASTQ join padding, default NNNNNNNN.
STR_OPT(join_padgapq)			// Q scores to use for FASTQ join padding, default IIIIIIII.
STR_OPT(sort)					// Sort order (length, size or other)
STR_OPT(uparseout)				// Output file, UPARSE tabbed format.
STR_OPT(sortedby)				// Sort order (length, size, other)
STR_OPT(padq)					//
STR_OPT(constax_report)			//
STR_OPT(sample)					//
STR_OPT(boot_subset)			//
STR_OPT(bitvec)				//
STR_OPT(sample_delim)			//
STR_OPT(mapout)			//
STR_OPT(rank)			//
STR_OPT(dbcutout)			//
STR_OPT(trimout)		//

UNS_OPT(stepwords,			8,			0,			UINT_MAX)	// -
UNS_OPT(bump,				50,			0,			100)		// -
UNS_OPT(rowlen,				80,			8,			UINT_MAX)	// Length of alignment row for -alnout.
UNS_OPT(mincodons,			20,			1,			UINT_MAX)	// Min number of codons in an ORF.
UNS_OPT(orfstyle,			(1+4),		0,			UINT_MAX)	// Options for ORF detection.
UNS_OPT(randseed,			1,			0,			UINT_MAX)	// Integer seed for random number generator.
UNS_OPT(minsize,			0,			0,			UINT_MAX)	// Min cluster size.
UNS_OPT(minuniquesize,		0,			0,			UINT_MAX)	// Min cluster size.
UNS_OPT(threads,			0,			0,			UINT_MAX)	// Number of threads (default: one per core).
UNS_OPT(maxhits,			0,			0,			UINT_MAX)	// Max number of hits to report for one query.
UNS_OPT(big,				100000,		0,			UINT_MAX)	// -
UNS_OPT(minseqlength,		8,			1,			UINT_MAX)	// Min sequence length (sortbysize output).
UNS_OPT(wordlength,			0,			1,			UINT_MAX)	// Word length for indexing.
UNS_OPT(dbstep,				0,			1,			UINT_MAX)	// Database step parameter.
UNS_OPT(slots,				0,			1,			UINT_MAX)	// Number of index slots.
UNS_OPT(dbaccelpct,			100,		1,			100)		// Database accel parameter.
UNS_OPT(posbits,			11,			1,			100)		// -
UNS_OPT(mosaic_minseg,		16,			1,			UINT_MAX)	// Min segment size.
UNS_OPT(maxdqm,				UINT_MAX,	0,			UINT_MAX)	// -
UNS_OPT(mindqt,				1,			1,			UINT_MAX)	// -

UNS_OPT(min_gene_length,	1200,		1,			UINT_MAX)	// -
UNS_OPT(max_gene_length,	2000,		1,			UINT_MAX)	// -

UNS_OPT(fastq_truncqual,	UINT_MAX,	0,			UINT_MAX)	// Truncate FASTQ read at first base with Q<=fastq_truncqual.
UNS_OPT(fastq_minqual,		UINT_MAX,	0,			UINT_MAX)	// Discard FASTQ read if any Q < minqual
UNS_OPT(fastq_minlen,		0,			0,			UINT_MAX)	// Min length of FASTQ read after filtering.
UNS_OPT(fastq_minovlen,		16,			4,			UINT_MAX)	// Min overlap of read pair for merging.
UNS_OPT(fastq_trunclen,		0,			0,			UINT_MAX)	// Truncate FASTQ at this length.
UNS_OPT(fastq_maxdiffs,		5,			0,			UINT_MAX)	// Max number of subs allowed in overlap for merging.
UNS_OPT(fastq_pctid,		90,			0,			UINT_MAX)	// Minimum alignment identity for merging
UNS_OPT(fastq_ascii,		33,			0,			UINT_MAX)	// ASCII base value for Q scores.
UNS_OPT(fastq_qmin,			0,			0,			UINT_MAX)	// Min Q score allowed by FASTQ format.
UNS_OPT(fastq_qmax,			42,			0,			UINT_MAX)	// Max Q score allowed by FASTQ format.
UNS_OPT(fastq_qmaxout,		42,			0,			UINT_MAX)	// Max Q score allowed by FASTQ format (output files).
UNS_OPT(fastq_tail,			4,			1,			UINT_MAX)	// Min number of bases in a FASTQ tail.
UNS_OPT(fastq_trunctail,	2,			0,			UINT_MAX)	// :doc Truncate tail
UNS_OPT(fastq_minmergelen,	0,			1,			UINT_MAX)	// Min length allowed for merged read.
UNS_OPT(fastq_maxmergelen,	0,			1,			UINT_MAX)	// Max length allowed for merged read.
UNS_OPT(fastq_stripleft,	0,			0,			UINT_MAX)	// Discard this number of bases at start of FASTQ read.
UNS_OPT(fastq_stripright,	0,			0,			UINT_MAX)	// Discard this number of bases at end of FASTQ read.
UNS_OPT(fastq_maxns,		UINT_MAX,	1,			UINT_MAX)	// :document
UNS_OPT(trunclen,			0,			1,			UINT_MAX)	// Truncate length.
UNS_OPT(padlen,				0,			1,			UINT_MAX)	// Pad with Ns if shorter.
UNS_OPT(stripleft,			0,			0,			UINT_MAX)	// Nr. positions to delete at start of sequence.
UNS_OPT(stripright,			0,			0,			UINT_MAX)	// Nr. positions to delete at end of sequence.
UNS_OPT(self_words_drop,	4,			1,			UINT_MAX)	// :document
UNS_OPT(flank,				8,			1,			UINT_MAX)	// :document
//UNS_OPT(minampsize,			8,			1,			UINT_MAX)	// :document
UNS_OPT(mincount,			UINT_MAX,	0,			UINT_MAX)	// :document

// uchime
UNS_OPT(mindiffs,			3,			1,			UINT_MAX)	// Mininum diffs to allow on each side of cross-over.
UNS_OPT(chunks,				4,			2,			UINT_MAX)	// Number of chunks.
UNS_OPT(minchunk,			64,			2,			UINT_MAX)	// Min chunk size.

UNS_OPT(maxaccepts,			1,			0,			UINT_MAX)	// Max number of accepts.
UNS_OPT(maxrejects,			32,			0,			UINT_MAX)	// Max number of rejects.

UNS_OPT(band,				16,			0,			UINT_MAX)	// Min band width.
UNS_OPT(hspw,				0,			1,			UINT_MAX)	// Word length for HSP finding.

UNS_OPT(minhsp,				16,			1,			UINT_MAX)	// Min length of HSP.
UNS_OPT(idprefix,			0,			0,			UINT_MAX)	// Min number of identical bases at start of alignment.
UNS_OPT(idsuffix,			0,			0,			UINT_MAX)	// Min number of identical bases at end of alignment.
UNS_OPT(maxdiffs,			UINT_MAX,	0,			UINT_MAX)	// Min number of diffs (subs+gaps) in the alignment.
UNS_OPT(maxsubs,			UINT_MAX,	0,			UINT_MAX)	// Min number of substitutions in the alignment.
UNS_OPT(maxgaps,			UINT_MAX,	0,			UINT_MAX)	// Min number of gapped columns in the alignment.
UNS_OPT(mincols,			0,			0,			UINT_MAX)	// Min number of columns in the alignment.
UNS_OPT(mindiffsa,			UINT_MAX,	0,			UINT_MAX)	// 
UNS_OPT(maxdiffsa,			UINT_MAX,	0,			UINT_MAX)	// 

UNS_OPT(maxseqlength,		50000,		1,			UINT_MAX)	// Max seqlence length (sortbylength).
UNS_OPT(topn,				UINT_MAX,	1,			UINT_MAX)	// Output only the first topn sequences (sortbysize, sortbylength).
UNS_OPT(uparse_maxhot,		32,			1,			UINT_MAX)	// Max nr of candidate parents to consider.
UNS_OPT(uparse_maxdrop,		8,			0,			UINT_MAX)	// -
UNS_OPT(uparse_maxdball,	100,		0,			UINT_MAX)	// -

UNS_OPT(boots,				100,		1,			UINT_MAX)	// -
UNS_OPT(default_size,		UINT_MAX,	1,			UINT_MAX)	// Default size for sortbysize.
UNS_OPT(long_target,		50000,		1,			UINT_MAX)	// -

UNS_OPT(fasta_cols,			80,			0,			UINT_MAX)	// Max length of a line when writing FASTA sequence, zero_array=all on one line.




UNS_OPT(maxstartdiffs,		4,			0,			UINT_MAX)	// -
UNS_OPT(maxenddiffs,		4,			0,			UINT_MAX)	// -




FLT_OPT(weak_id,			0.0,		0.0,		1.0)		// Min identity for weak hit.
FLT_OPT(weak_evalue,		10.0,		0.0,		DBL_MAX)	// Max E-value for weak hit.
FLT_OPT(ka_gapped_lambda,	0.0,		0.0,		DBL_MAX)	// Lambda for gapped alignments.
FLT_OPT(ka_ungapped_lambda,	0.0,		0.0,		DBL_MAX)	// Lambda for ungapped alignments.
FLT_OPT(ka_gapped_k,		0.0,		0.0,		DBL_MAX)	// K for gapped alignments.
FLT_OPT(ka_ungapped_k,		0.0,		0.0,		DBL_MAX)	// K for ungapped alignments.
FLT_OPT(ka_dbsize,			3e9,		1.0,		DBL_MAX)	// Effective database size.
FLT_OPT(accel,				0.8,		0.0,		1.0)		// Accel parameter for UBLAST.

FLT_OPT(uparse_match,		0,			-1000.0,	1000.0)		// Match score for UPARSE-REF d.p.
FLT_OPT(uparse_mismatch,	-1.0,		-1000.0,	1000.0)		// Mismatch score for UPARSE-REF d.p.
FLT_OPT(uparse_break,		-3.0,		-1000.0,	1000.0)		// Chimeric breakpoint score for UPARSE-REF d.p.

FLT_OPT(id,					0.0,		0.0,		1.0)		// Min identity (0.0 to 1.0) for search.
FLT_OPT(mid,				0.0,		0.0,		1.0)		// Min match identity (0.0 to 1.0) for search.
FLT_OPT(evalue,				10.0,		0.0,		DBL_MAX)	// Max E-value.
FLT_OPT(minqt,				0.0,		0.0,		DBL_MAX)	// Min query/target length ratio.
FLT_OPT(maxqt,				DBL_MAX,	0.0,		DBL_MAX)	// Max query/target length ratio.
FLT_OPT(minsl,				0.0,		0.0,		DBL_MAX)	// Min shorter/longer length ratio.
FLT_OPT(maxsl,				DBL_MAX,	0.0,		DBL_MAX)	// Max shorter/longer length ratio.
FLT_OPT(maxid,				1.0,		0.0,		1.0)		// Max identity (0.0 to 1.0) for hit.
FLT_OPT(min_sizeratio,		0.0,		0.0,		DBL_MAX)	// Min size ratio.

FLT_OPT(match,				1.0,		0.0,		DBL_MAX)	// Match score.
FLT_OPT(mismatch,			-2.0,		0.0,		DBL_MAX)	// Mismatch score.
FLT_OPT(xdrop_u,			16.0,		0.0,		DBL_MAX)	// X-drop parameter for ungapped alignments.
FLT_OPT(xdrop_g,			32.0,		0.0,		DBL_MAX)	// X-drop parameter for gapped alignments.
FLT_OPT(xdrop_nw,			8.0,		0.0,		DBL_MAX)	// X-drop parameter for HSP-finding.

FLT_OPT(minh,				0.35,		0.0,		DBL_MAX)	// Min chimera score for yes
FLT_OPT(xn,					8.0,		0.0,		DBL_MAX)	// No vote wight.
FLT_OPT(dn,					1.4,		0.0,		DBL_MAX)	// Pseudo-count prior.
FLT_OPT(xa,					1.0,		0.0,		DBL_MAX)	// Abstain vote weight.
FLT_OPT(mindiv,				1.0,		0.0,		100.0)		// Min chimera divergence.

FLT_OPT(uparse_annot_maxdivqm,1.0,		0.0,		100.0)		// Max divergence between query and model.

FLT_OPT(fastq_maxee,		DBL_MAX,	0.0,		DBL_MAX)	// Max expected errors in FASTQ read.
FLT_OPT(fastq_maxee_rate,	DBL_MAX,	0.0,		1.0)		// Max expected errors per base in FASTQ read.

FLT_OPT(abskew,				2,			0.0,		100.0)		// Min abundance skew.

FLT_OPT(query_cov,			0.0,		0.0,		1.0)		// Min fraction of query covered by alignment.
FLT_OPT(target_cov,			0.0,		0.0,		1.0)		// Min fraction of target covered by alignment.
FLT_OPT(max_query_cov,		0.0,		0.0,		1.0)		// Min fraction of query covered by alignment.
FLT_OPT(max_target_cov,		0.0,		0.0,		1.0)		// Min fraction of target covered by alignment.

FLT_OPT(lopen,				10.0,		0.0,		DBL_MAX)	// Gap open penalty for local alignments.
FLT_OPT(lext,				1.0,		0.0,		DBL_MAX)	// Gap extend penalty for local alignments.

FLT_OPT(termid,				0.0,		0.0,		1.0)		// Terminate search if identity < termid.
FLT_OPT(termidd,			0.0,		0.0,		1.0)		// Terminate search if (maxid-id) > termidd.

FLT_OPT(orient_wordx,		8.0,		0.0,		DBL_MAX)	// - 
FLT_OPT(orient_strandx,		4.0,		0.0,		DBL_MAX)	// - 
FLT_OPT(sintax_cutoff,		0.8,		0.0,		1.0)		// - 
FLT_OPT(unoise_alpha,		2.0,		1.0,		999)



FLT_OPT(maj,				0.51,		0.0,		1.0)



FLAG_OPT(quiet)				// Show only errors and warning messages on the terminal.
FLAG_OPT(quicksort)			// -
FLAG_OPT(kmerid)			// -
FLAG_OPT(self)				// Ignore hits with same label.
FLAG_OPT(notself)			// Ignore hits with different labels.
FLAG_OPT(selfid)			// Ignore hits to identical sequence.
FLAG_OPT(fulldp)			// 
FLAG_OPT(logmemgrows)		// -
FLAG_OPT(trunclabels)		// Truncate sequence labels.
FLAG_OPT(notrunclabels)		// Truncate sequence labels.
FLAG_OPT(despacelabels)		// Replace white space in labels with _
FLAG_OPT(verbose)			// -
FLAG_OPT(log_query)			// -
FLAG_OPT(sizein)			// Count "size=N;" annotations towards cluster sizes.
FLAG_OPT(sizeout)			// Add "size=N;" annotations to cluster centroid labels.
FLAG_OPT(hardmask)			// Mask by replacing letters with wildcard (N or X).
FLAG_OPT(validate)			// -
FLAG_OPT(uc_hitsonly)		// Omit N records in uc file
FLAG_OPT(gaforce)			// Force global alignment even if HSPs are below threshold.
FLAG_OPT(random_top_hit)	// Report only top hit (ties broken randomly).
FLAG_OPT(top_hit_only)		// Report only top hit (ties broken systematically).
FLAG_OPT(bottom_hit_only)	// Report only bottom hit (ties broken systematically).
FLAG_OPT(top_hits_only)		// Report only top hits (including ties).

FLAG_OPT(fastq_nostagger)			// Don't allow staggered merges

FLAG_OPT(leftjust)			// Reject hit if terminal gaps at start of alignment.
FLAG_OPT(rightjust)			// Reject hit if terminal gaps at end of alignment.
FLAG_OPT(notermgapsq)		// Reject hit if terminal gaps in query

FLAG_OPT(log_u)					// -
FLAG_OPT(log_sortedu)			// -
FLAG_OPT(log_searcher_alns)		// -
FLAG_OPT(output_no_hits)		// Report queries with no hits.
FLAG_OPT(fastapairs_dots)		// Replace identities with dots in query sequence.
FLAG_OPT(show_termgaps)			// Show terminal gaps in alnout.
FLAG_OPT(end_of_row)			// -
FLAG_OPT(cover_query)			// Report sufficient hits to cover query sequence.
FLAG_OPT(mosaic)				// Report mosaic hits.
FLAG_OPT(strict_newick)			// Re-format labels to confirm with Newick specification.
FLAG_OPT(fastq_logvaln)			// Verbose mergepairs vertical alignment (log file).
FLAG_OPT(fastq_eeout)			// 
FLAG_OPT(constax)				// 
FLAG_OPT(keepgaps)				// 
FLAG_OPT(interleaved)			//
FLAG_OPT(orf_plusonly)			//
//FLAG_OPT(uchime_sensitive)		//
FLAG_OPT(merge_annot)			//
FLAG_OPT(merge_ambig)			//
FLAG_OPT(offby1)			//
FLAG_OPT(allxch)			//
FLAG_OPT(maxskew)			//
FLAG_OPT(bimeras_only)		//
FLAG_OPT(overlap_only)		//
FLAG_OPT(fastq_noguess)		//
FLAG_OPT(fastq_forceq)		//
FLAG_OPT(tax_prod)		//
FLAG_OPT(ktop)		//
FLAG_OPT(ignore_label_mismatches)
FLAG_OPT(tov)
FLAG_OPT(fasta)
FLAG_OPT(fastq)
FLAG_OPT(force_unique_labels)
FLAG_OPT(allow_digits)

VECTOR_OPT(fastq_mergepairs)
VECTOR_OPT(fastq_filter)
VECTOR_OPT(reverse)

#undef FLAG_OPT
#undef UNS_OPT
#undef FLT_OPT
#undef STR_OPT
#undef VECTOR_OPT
