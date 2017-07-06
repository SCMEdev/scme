
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H0 X
# H0 X   libAtoms+QUIP: atomistic simulation library
# H0 X
# H0 X   Portions of this code were written by
# H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# H0 X
# H0 X   Copyright 2006-2010.
# H0 X
# H0 X   These portions of the source code are released under the GNU General
# H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# H0 X
# H0 X   If you would like to license the source code under different terms,
# H0 X   please contact Gabor Csanyi, gabor@csanyi.net
# H0 X
# H0 X   Portions of this code were written by Noam Bernstein as part of
# H0 X   his employment for the U.S. Government, and are not subject
# H0 X   to copyright in the USA.
# H0 X
# H0 X
# H0 X   When using this software, please cite the following reference:
# H0 X
# H0 X   http://www.libatoms.org
# H0 X
# H0 X  Additional contributions by
# H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# H0 X
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

.PHONY: default all clean default_scme

include Makefile.${QUIP_ARCH}
include Makefile.inc
include Makefile.rules

THIRDPARTY_LIBS := libthirdparty.a
ifeq (${HAVE_FX},1)
THIRDPARTY_LIBS += libfx.a
endif

ifeq (${HAVE_SCME},1)
THIRDPARTY_LIBS += libscme.a
endif

ifeq (${HAVE_MTP},1)
THIRDPARTY_LIBS += libmtp.a
endif

default: ${THIRDPARTY_LIBS}
all: default

LIBTHIRDPARTY_F95_FILES = \
	cube_tools

LIBTHIRDPARTY_F77_FILES= 

LIBTHIRDPARTY_C_FILES =


LIBTHIRDPARTY_F77_SOURCES = ${addsuffix .f, ${LIBTHIRDPARTY_F77_FILES}}
LIBTHIRDPARTY_F77_OBJS = ${addsuffix .o, ${LIBTHIRDPARTY_F77_FILES}}
LIBTHIRDPARTY_F95_SOURCES = ${addsuffix .f95, ${LIBTHIRDPARTY_F95_FILES}}
LIBTHIRDPARTY_F95_OBJS = ${addsuffix .o, ${LIBTHIRDPARTY_F95_FILES}}
LIBTHIRDPARTY_C_SOURCES = ${addsuffix .c, ${LIBTHIRDPARTY_C_FILES}}
LIBTHIRDPARTY_C_OBJS = ${addsuffix .o, ${LIBTHIRDPARTY_C_FILES}}

LIBTHIRDPARTY_SOURCES = ${LIBTHIRDPARTY_F77_SOURCES} ${LIBTHIRDPARTY_F95_SOURCES} ${LIBTHIRDPARTY_C_SOURCES}
LIBTHIRDPARTY_OBJS = ${LIBTHIRDPARTY_F77_OBJS} ${LIBTHIRDPARTY_F95_OBJS} ${LIBTHIRDPARTY_C_OBJS}

libthirdparty.a: ${LIBTHIRDPARTY_OBJS}
ifneq (${LIBTOOL},)
	${LIBTOOL} -o libthirdparty.a ${LIBTHIRDPARTY_OBJS}
else
	${AR} ${AR_ADD} libthirdparty.a $?
endif

libfx.a: nttm3f
	cp ${PWD}/src/ThirdParty/nttm3f/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/ThirdParty/nttm3f -I${PWD} -I${PWD}/arch $@

libscme.a: default_scme
	cp ${PWD}/src/ThirdParty/scme/quip/Makefile.SCME ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/ThirdParty/scme/src -I${PWD} -I${PWD}/arch $@

libmtp.a: mtp
	cp ${PWD}/src/ThirdParty/mtp/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/ThirdParty/mtp -I${PWD} -I${PWD}/arch $@

clean:
	rm -f libthirdparty.a libfx.a libscme.a libmtp.a *.o *.mod *.mod.save *.mod.txt *.mod.save.txt *.fpp *.f90doc

allclean: clean 

install:

${LIBTHIRDPARTY_C_OBJS}  : %.o : %.c
	${CC} ${CFLAGS} -c $< 


