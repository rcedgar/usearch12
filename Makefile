BINDIR := ../bin
OBJDIR := o
BINPATH := $(BINDIR)/usearch12

CPPFLAGS := $(CPPFLAGS) -DNDEBUG -pthread

CC = ccache gcc
CFLAGS := $(CFLAGS) -O3 -fopenmp -ffast-math -msse -mfpmath=sse

CXX = ccache g++
CXXFLAGS := $(CXXFLAGS) -O3 -fopenmp -ffast-math -msse -mfpmath=sse

UNAME_S := $(shell uname -s)
LDFLAGS := $(LDFLAGS) -O3 -fopenmp -pthread -lpthread
ifeq ($(UNAME_S),Linux)
    LDFLAGS += -static
endif

HDRS = \
  accepter.h \
  aligner.h \
  alignresult.h \
  alnheuristics.h \
  alnparams.h \
  alpha.h \
  alphainfo.h \
  bitmapsearcher.h \
  bitvec.h \
  bitvec64.h \
  chainer.h \
  chainer1.h \
  chimehit.h \
  chunksearcher.h \
  closedrefsink.h \
  clustersink.h \
  cmd.h \
  cmds.h \
  constaxf.h \
  constaxsink.h \
  constaxstr.h \
  count.h \
  countsort.h \
  cpplock.h \
  crc32.h \
  dbhitsink.h \
  dbtype.h \
  deflate.h \
  deparser.h \
  derep.h \
  derepresult.h \
  diagbox.h \
  duster.h \
  estats.h \
  evalue.h \
  fastaseqsource.h \
  fastq.h \
  fastqseqsource.h \
  fileseqsource.h \
  filetype.h \
  finger.h \
  fragaligner.h \
  genefinder.h \
  globalaligner.h \
  gobuff.h \
  gzguts.h \
  hitmgr.h \
  hitsink.h \
  hsp.h \
  hspfinder.h \
  inffast.h \
  inffixed.h \
  inflate.h \
  inftrees.h \
  label.h \
  linereader.h \
  localaligner.h \
  localaligner2.h \
  lockobj.h \
  lockobjs.h \
  mask.h \
  merge.h \
  mergeglobals.h \
  mx.h \
  myopts.h \
  myutils.h \
  obj.h \
  objmgr.h \
  objtype.h \
  objtypes.h \
  orffinder.h \
  otutab.h \
  otutabsink.h \
  outputsink.h \
  pathinfo.h \
  pcb.h \
  primes.h \
  quarts.h \
  searcher.h \
  segmask.h \
  seqdb.h \
  seqdbsearcher.h \
  seqdbseqsource.h \
  seqhash.h \
  seqinfo.h \
  seqsource.h \
  sintaxsearcher.h \
  sort.h \
  sparsemx.h \
  strdict.h \
  tax.h \
  taxy.h \
  terminator.h \
  tracebit.h \
  tree.h \
  trees.h \
  udbdata.h \
  udbfile.h \
  udbparams.h \
  udbsearcher.h \
  udbusortedsearcher.h \
  ufenum.h \
  uparsesink.h \
  upclustersink.h \
  userfields.h \
  wordcounter.h \
  xdpmem.h \
  xtype.h \
  zconf.h \
  zlib.h \
  zutil.h \

OBJS = \
  $(OBJDIR)/accepter.o \
  $(OBJDIR)/alignresult.o \
  $(OBJDIR)/alnheuristics.o \
  $(OBJDIR)/alnout.o \
  $(OBJDIR)/alnparams.o \
  $(OBJDIR)/alpha.o \
  $(OBJDIR)/alpha2.o \
  $(OBJDIR)/alphainfo.o \
  $(OBJDIR)/arscorer.o \
  $(OBJDIR)/bimeradp.o \
  $(OBJDIR)/bitmapsearcher.o \
  $(OBJDIR)/bitvec.o \
  $(OBJDIR)/bitvec64.o \
  $(OBJDIR)/blast6out.o \
  $(OBJDIR)/blosum62.o \
  $(OBJDIR)/chainer1.o \
  $(OBJDIR)/chainer.o \
  $(OBJDIR)/chimehit.o \
  $(OBJDIR)/chunksearcher.o \
  $(OBJDIR)/closedrefsink.o \
  $(OBJDIR)/clusterfast.o \
  $(OBJDIR)/clustermt.o \
  $(OBJDIR)/clustersink.o \
  $(OBJDIR)/clustersmallmem.o \
  $(OBJDIR)/cmd.o \
  $(OBJDIR)/cmdline.o \
  $(OBJDIR)/comppath.o \
  $(OBJDIR)/constaxf.o \
  $(OBJDIR)/constaxsink.o \
  $(OBJDIR)/constaxstr.o \
  $(OBJDIR)/countsort.o \
  $(OBJDIR)/dbhitsink.o \
  $(OBJDIR)/deparser.o \
  $(OBJDIR)/derepfull.o \
  $(OBJDIR)/derepresult.o \
  $(OBJDIR)/diagbox.o \
  $(OBJDIR)/dustmask.o \
  $(OBJDIR)/estats.o \
  $(OBJDIR)/evalue.o \
  $(OBJDIR)/fastaseqsource.o \
  $(OBJDIR)/fastqfilter2.o \
  $(OBJDIR)/fastxgetsamplenames.o \
  $(OBJDIR)/fastmask.o \
  $(OBJDIR)/fastq.o \
  $(OBJDIR)/fastqfilter.o \
  $(OBJDIR)/fastqjoin.o \
  $(OBJDIR)/fastqmerge.o \
  $(OBJDIR)/fastqseqsource.o \
  $(OBJDIR)/fastxtruncate.o \
  $(OBJDIR)/filetype.o \
  $(OBJDIR)/fileseqsource.o \
  $(OBJDIR)/findgene.o \
  $(OBJDIR)/finger.o \
  $(OBJDIR)/genefinder.o \
  $(OBJDIR)/getcmd.o \
  $(OBJDIR)/getfastqs.o \
  $(OBJDIR)/getglobalhsps.o \
  $(OBJDIR)/gethsps.o \
  $(OBJDIR)/getuniquelettercount.o \
  $(OBJDIR)/globalaligner.o \
  $(OBJDIR)/globalalignmem.o \
  $(OBJDIR)/makeclustersearcher.o \
  $(OBJDIR)/mergealign.o \
  $(OBJDIR)/mergelogvaln.o \
  $(OBJDIR)/gzipfileio.o \
  $(OBJDIR)/hitmgr.o \
  $(OBJDIR)/hspfinder.o \
  $(OBJDIR)/interp.o \
  $(OBJDIR)/json.o \
  $(OBJDIR)/label.o \
  $(OBJDIR)/linereader.o \
  $(OBJDIR)/lnfrac.o \
  $(OBJDIR)/loaddb.o \
  $(OBJDIR)/localaligner.o \
  $(OBJDIR)/localaligner2.o \
  $(OBJDIR)/localmulti.o \
  $(OBJDIR)/lockobj.o \
  $(OBJDIR)/logaln.o \
  $(OBJDIR)/make3way.o \
  $(OBJDIR)/makedbsearcher.o \
  $(OBJDIR)/makeudb.o \
  $(OBJDIR)/mask.o \
  $(OBJDIR)/mergepair.o \
  $(OBJDIR)/mergepost.o \
  $(OBJDIR)/mergepre.o \
  $(OBJDIR)/mergestats.o \
  $(OBJDIR)/mergethread.o \
  $(OBJDIR)/mx.o \
  $(OBJDIR)/myutils.o \
  $(OBJDIR)/newick.o \
  $(OBJDIR)/objmgr.o \
  $(OBJDIR)/fragaligner.o \
  $(OBJDIR)/orffinder.o \
  $(OBJDIR)/orient.o \
  $(OBJDIR)/otutab.o \
  $(OBJDIR)/otutabsink.o \
  $(OBJDIR)/outputsink.o \
  $(OBJDIR)/outputuc.o \
  $(OBJDIR)/pathinfo.o \
  $(OBJDIR)/pattern.o \
  $(OBJDIR)/pcb.o \
  $(OBJDIR)/prime.o \
  $(OBJDIR)/quarts.o \
  $(OBJDIR)/search.o \
  $(OBJDIR)/searchcmd.o \
  $(OBJDIR)/searcher.o \
  $(OBJDIR)/segmaskseq.o \
  $(OBJDIR)/seqdb.o \
  $(OBJDIR)/seqdbfromfasta.o \
  $(OBJDIR)/seqdbio.o \
  $(OBJDIR)/seqdbsearcher.o \
  $(OBJDIR)/seqdbseqsource.o \
  $(OBJDIR)/seqhash.o \
  $(OBJDIR)/seqinfo.o \
  $(OBJDIR)/seqsource.o \
  $(OBJDIR)/setnucmx.o \
  $(OBJDIR)/sintaxsummary.o \
  $(OBJDIR)/sintaxsearcher.o \
  $(OBJDIR)/staralign.o \
  $(OBJDIR)/strdict.o \
  $(OBJDIR)/substmx.o \
  $(OBJDIR)/tax.o \
  $(OBJDIR)/taxy.o \
  $(OBJDIR)/terminator.o \
  $(OBJDIR)/tracebackbitmem.o \
  $(OBJDIR)/tree.o \
  $(OBJDIR)/treefromagg.o \
  $(OBJDIR)/treetofile.o \
  $(OBJDIR)/uchime3denovo.o \
  $(OBJDIR)/udbbuild.o \
  $(OBJDIR)/udb2bitvec.o \
  $(OBJDIR)/udbusortedsearcherbig.o \
  $(OBJDIR)/udbdata.o \
  $(OBJDIR)/udbio.o \
  $(OBJDIR)/udbparams.o \
  $(OBJDIR)/udbsearcher.o \
  $(OBJDIR)/udbusortedsearcher.o \
  $(OBJDIR)/ungappedblast.o \
  $(OBJDIR)/unoise3.o \
  $(OBJDIR)/uparsesink.o \
  $(OBJDIR)/uparsedp.o \
  $(OBJDIR)/uparsepretty.o \
  $(OBJDIR)/upclustersink.o \
  $(OBJDIR)/usearch_main.o \
  $(OBJDIR)/userout.o \
  $(OBJDIR)/viterbifastbandmem.o \
  $(OBJDIR)/viterbifastmem.o \
  $(OBJDIR)/wordcounter.o \
  $(OBJDIR)/wordparams.o \
  $(OBJDIR)/xdropalignmem.o \
  $(OBJDIR)/xdropbwdmem.o \
  $(OBJDIR)/xdropbwdsplit.o \
  $(OBJDIR)/xdropfwdmem.o \
  $(OBJDIR)/xdropfwdsplit.o \
  $(OBJDIR)/adler32.o \
  $(OBJDIR)/crc32.o \
  $(OBJDIR)/deflate.o \
  $(OBJDIR)/gzlib.o \
  $(OBJDIR)/gzread.o \
  $(OBJDIR)/infback.o \
  $(OBJDIR)/inffast.o \
  $(OBJDIR)/inflate.o \
  $(OBJDIR)/inftrees.o \
  $(OBJDIR)/trees.o \
  $(OBJDIR)/zutil.o \

.PHONY: clean

$(BINPATH) : $(BINDIR)/ $(OBJDIR)/ $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(BINPATH)
	strip -d $(BINPATH)

$(OBJDIR)/ :
	mkdir -p $(OBJDIR)/

$(BINDIR)/ :
	mkdir -p $(BINDIR)/

$(OBJDIR)/%.o : %.c $(HDRS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(OBJDIR)/%.o : %.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR)/ $(BINPATH)
