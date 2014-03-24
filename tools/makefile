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
BAMTOOLSDIR = ../external/BamTools
BAMTOOLSLIB = libbamtools.a

INCLUDE = -I. -I../include -I$(BAMDIR) -I$(ZDIR) -I$(BAMTOOLSDIR)
CC = /usr/bin/g++
CXX = /usr/bin/g++
CXXFLAGS = -O3 -DDEBUG_CHECKS $(INCLUDE)
LDFLAGS = -O3 -L$(BAMDIR) -lbam -L$(ZDIR) -lz -L$(BAMTOOLSDIR) -lbamtools

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

$(BAMTOOLSDIR)/$(BAMTOOLSLIB) :
		$(MAKE) --no-print-directory -C $(BAMTOOLSDIR) $(BAMTOOLSLIB)

#
# Build tools
#

clustermatepairs : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./MatePairEM.o\
		./asa136.o\
		./asa241.o\
		./clustermatepairs.o\
		$(BAMTOOLSDIR)/$(BAMTOOLSLIB)\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

setcover: \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./setcover.o\
		$(BAMTOOLSDIR)/$(BAMTOOLSLIB)\
		$(ZDIR)/$(ZLIB)

initsplitalign: \
		./AlignmentStream.o\
		./FastaIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./ReadIndex.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./initsplitalign.o\
		$(BAMTOOLSDIR)/$(BAMTOOLSLIB)\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

dosplitalign: \
		./AlignmentStream.o\
		./FastaIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./ReadIndex.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./ReadStream.o\
		./dosplitalign.o\
		$(BAMTOOLSDIR)/$(BAMTOOLSLIB)\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

evalsplitalign: \
		./AlignmentIndex.o\
		./FastaIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./ReadIndex.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./ReadStream.o\
		./evalsplitalign.o\
		$(BAMTOOLSDIR)/$(BAMTOOLSLIB)\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

splitseq: \
		./AlignmentStream.o\
		./FastaIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./ReadIndex.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./PathFinder.o\
		./splitseq.o\
		$(BAMTOOLSDIR)/$(BAMTOOLSLIB)\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

calccov: \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./calccov.o\
		$(BAMTOOLSDIR)/$(BAMTOOLSLIB)\
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

estislands : \
		./estislands.o\
		./EstCatalog.o\
		$(BAMTOOLSDIR)/$(BAMTOOLSLIB)\
		$(ZDIR)/$(ZLIB)

localalign : \
		./localalign.o\
		./SimpleAligner.o\
		./DebugCheck.o\
		./Common.o

matealign : \
		./matealign.o\
		./Common.o\
		./DebugCheck.o\
		./ReadStream.o\
		./Sequences.o\
		./SimpleAligner.o

tools : clustermatepairs setcover calccov estislands initsplitalign dosplitalign evalsplitalign localalign splitseq matealign
	

clean : 
		rm -f *.o
		$(MAKE) --no-print-directory -C $(ZDIR) clean
		$(MAKE) --no-print-directory -C $(BAMDIR) clean
