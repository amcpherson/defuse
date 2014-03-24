###################################################
#
# Makefile for tools
# Creator [Xcode -> Makefile Ver: Feb 14 2007 09:18:41]
# Created: [Tue Aug 10 13:19:17 2010]
#
###################################################

#
# Macros
#

BAMDIR = ../external/samtools-0.1.8
BAMLIB = libbam.a
ZDIR = ../external/zlib-1.2.5
ZLIB = libz.a

INCLUDE = -I. -I../include -I$(BAMDIR) -I$(ZDIR)
CC = /usr/bin/g++
CXX = /usr/bin/g++
CXXFLAGS = -O2 -DDEBUG_CHECKS $(INCLUDE)
LDFLAGS = -O2 -L$(BAMDIR) -lbam -L$(ZDIR) -lz

all: tools $(BAMDIR)/samtools

#
# Build samtools
#

$(BAMDIR)/samtools :
		$(MAKE) --no-print-directory -C $(BAMDIR) samtools

# 
# Build libs
# 

$(ZDIR)/$(ZLIB) : 
		$(MAKE) --no-print-directory -C $(ZDIR) $(ZLIB)

$(BAMDIR)/$(BAMLIB) : 
		$(MAKE) --no-print-directory -C $(BAMDIR) $(BAMLIB)

#
# Build tools
#

connectivity : \
		./Common.o\
		./QuaterpartiteEdges.o\
		./DebugCheck.o\
		./Parsers.o\
		./connectivity.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

remapsplices: \
		./Common.o\
		./CompactBreakRegion.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./remapsplices.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

overlapclusters : \
		./Common.o\
		./CompactBreakRegion.o\
		./DebugCheck.o\
		./overlapclusters.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

maximalvalidclusters : \
		./Clusters.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./MatePairCombinatorial.o\
		./maximalvalidclusters.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

setcover: \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./ReweightSetCover.o\
		./setcover.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

splitseq: \
		./AlignmentIndex.o\
		./ReadIndex.o\
		./FastaIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./PathFinder.o\
		./SplitFusion.o\
		./splitseq.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

denovoseq: \
		./AlignmentIndex.o\
		./ReadIndex.o\
		./FastaIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./PathFinder.o\
		./DeNovoFusion.o\
		./denovoseq.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

calccov: \
		./AlignmentIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./calccov.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

estislands : \
		./estislands.o\
		./EstCatalog.o\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

tools : connectivity remapsplices overlapclusters maximalvalidclusters setcover splitseq denovoseq calccov estislands
	

clean : 
		rm -f *.o
		$(MAKE) --no-print-directory -C $(ZDIR) clean
		$(MAKE) --no-print-directory -C $(BAMDIR) clean
