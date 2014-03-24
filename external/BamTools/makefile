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

INCLUDE = -I.
CXX = /usr/bin/g++
CXXFLAGS = -D_FILE_OFFSET_BITS=64 -DBAMTOOLS_API_LIBRARY -O3 -DNDEBUG -Wall -fPIC $(INCLUDE)

all: libbamtools.a

#
# Build lib
#

LOBJS=\
	api/BamAlignment.o \
	api/BamMultiReader.o \
	api/BamReader.o \
	api/BamWriter.o \
	api/SamHeader.o \
	api/SamProgram.o \
	api/SamProgramChain.o \
	api/SamReadGroup.o \
	api/SamReadGroupDictionary.o \
	api/SamSequence.o \
	api/SamSequenceDictionary.o \
	api/internal/BamDeviceFactory_p.o \
	api/internal/BamException_p.o \
	api/internal/BamFile_p.o \
	api/internal/BamFtp_p.o \
	api/internal/BamHeader_p.o \
	api/internal/BamHttp_p.o \
	api/internal/BamIndexFactory_p.o \
	api/internal/BamMultiReader_p.o \
	api/internal/BamPipe_p.o \
	api/internal/BamRandomAccessController_p.o \
	api/internal/BamReader_p.o \
	api/internal/BamStandardIndex_p.o \
	api/internal/BamToolsIndex_p.o \
	api/internal/BamWriter_p.o \
	api/internal/BgzfStream_p.o \
	api/internal/ILocalIODevice_p.o \
	api/internal/SamFormatParser_p.o \
	api/internal/SamFormatPrinter_p.o \
	api/internal/SamHeaderValidator_p.o

libbamtools.a:$(LOBJS)
	$(AR) cr $@ $(LOBJS)

clean : 
	rm -f $(LOBJS) libbamtools.a

