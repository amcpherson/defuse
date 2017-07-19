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
BAMTOOLSDIR = ../external/BamTools-2.0/
BAMTOOLSAPIDIR = ../external/BamTools-2.0/api

INCLUDE = -I. -I../include -I$(BAMDIR) -I$(ZDIR) -I$(BAMTOOLSDIR) -I$(BAMTOOLSAPIDIR)
CC = g++
CXX = g++
CXXFLAGS = -O3 -g -DDEBUG_CHECKS $(INCLUDE)
LDFLAGS = -O3 -g -L$(ZDIR) -lz

all: tools

# 
# Build libs
# 

$(ZDIR)/$(ZLIB) :
		$(MAKE) --no-print-directory -C $(ZDIR) $(ZLIB)

$(BAMDIR)/$(BAMLIB) :
		$(MAKE) --no-print-directory -C $(BAMDIR) $(BAMLIB)

#
# BamTools
#

BAMTOOLSSRC = \
		$(BAMTOOLSAPIDIR)/BamAlignment.o\
		$(BAMTOOLSAPIDIR)/BamMultiReader.o\
		$(BAMTOOLSAPIDIR)/BamReader.o\
		$(BAMTOOLSAPIDIR)/BamWriter.o\
		$(BAMTOOLSAPIDIR)/SamHeader.o\
		$(BAMTOOLSAPIDIR)/SamProgram.o\
		$(BAMTOOLSAPIDIR)/SamProgramChain.o\
		$(BAMTOOLSAPIDIR)/SamReadGroup.o\
		$(BAMTOOLSAPIDIR)/SamReadGroupDictionary.o\
		$(BAMTOOLSAPIDIR)/SamSequence.o\
		$(BAMTOOLSAPIDIR)/SamSequenceDictionary.o\
		$(BAMTOOLSAPIDIR)/internal/bam/BamHeader_p.o\
		$(BAMTOOLSAPIDIR)/internal/bam/BamMultiReader_p.o\
		$(BAMTOOLSAPIDIR)/internal/bam/BamRandomAccessController_p.o\
		$(BAMTOOLSAPIDIR)/internal/bam/BamReader_p.o\
		$(BAMTOOLSAPIDIR)/internal/bam/BamWriter_p.o\
		$(BAMTOOLSAPIDIR)/internal/index/BamIndexFactory_p.o\
		$(BAMTOOLSAPIDIR)/internal/index/BamStandardIndex_p.o\
		$(BAMTOOLSAPIDIR)/internal/index/BamToolsIndex_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/BamDeviceFactory_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/BamFile_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/BamFtp_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/BamHttp_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/BamPipe_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/BgzfStream_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/ByteArray_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/HostAddress_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/HostInfo_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/HttpHeader_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/ILocalIODevice_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/RollingBuffer_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/TcpSocket_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/TcpSocketEngine_p.o\
		$(BAMTOOLSAPIDIR)/internal/io/TcpSocketEngine_unix_p.o\
		$(BAMTOOLSAPIDIR)/internal/sam/SamFormatParser_p.o\
		$(BAMTOOLSAPIDIR)/internal/sam/SamFormatPrinter_p.o\
		$(BAMTOOLSAPIDIR)/internal/sam/SamHeaderValidator_p.o\
		$(BAMTOOLSAPIDIR)/internal/utils/BamException_p.o

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
		./clustermatepairs.o

setcover: \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./setcover.o

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
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

evalsplitalign: \
		./FastaIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./ReadIndex.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./ReadStream.o\
		./evalsplitalign.o\
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
		$(BAMDIR)/$(BAMLIB)\
		$(ZDIR)/$(ZLIB)

calccov: \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./calccov.o

estislands : \
		./estislands.o\
		./EstCatalog.o

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

bamfastq : \
		./Common.o\
		./DebugCheck.o\
		./bamfastq.o\
		$(BAMTOOLSSRC)\
		$(ZDIR)/$(ZLIB)

tools : $(ZDIR)/$(ZLIB) $(BAMDIR)/$(BAMLIB) clustermatepairs setcover calccov estislands dosplitalign evalsplitalign localalign splitseq matealign bamfastq
	

clean : 
		rm -f *.o
		$(MAKE) --no-print-directory -C $(ZDIR) clean
