# Makefile for zlib
# Copyright (C) 1995-2010 Jean-loup Gailly.
# For conditions of distribution and use, see copyright notice in zlib.h

# To compile and test, type:
#    ./configure; make test
# Normally configure builds both a static and a shared library.
# If you want to build just a static library, use: ./configure --static

# To use the asm code, type:
#    cp contrib/asm?86/match.S ./match.S
#    make LOC=-DASMV OBJA=match.o

# To install /usr/local/lib/libz.* and /usr/local/include/zlib.h, type:
#    make install
# To install in $HOME instead of /usr/local, use:
#    make install prefix=$HOME

CC = gcc
CFLAGS = -O3 -D_LARGEFILE64_SOURCE=1

AR=ar rc

OBJS = adler32.o compress.o crc32.o deflate.o gzclose.o gzlib.o gzread.o \
	gzwrite.o infback.o inffast.o inflate.o inftrees.o trees.o uncompr.o zutil.o

libz.a: $(OBJS)
	$(AR) $@ $(OBJS)

clean: 
	rm -f $(OBJS) libz.a
