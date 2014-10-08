# Request latest gcc!
# . /opt/rh/devtoolset-2/enable

VR       := 1
VECHOST  := -mavx -vec-report=${VR}
VECMIC   := -mmic -vec-report=${VR}

### To disable vectorization set USER_CXXFLAGS := -no-simd -no-vec
# Setting only one of the above has little effect.
# Note, this also screws-up prefetching so it's a lousy deal.

# -opt-prefetch-distance=64,8
OPT      := -O3

CPPFLAGS := -I. ${USER_CPPFLAGS} ${DEFS}
CXXFLAGS := ${OPT} -openmp -std=gnu++0x ${USER_CXXFLAGS}

LDFLAGS  := ${USER_LDFLAGS}

TGTS     := mkFit

EXES     := ${TGTS} $(addsuffix -mic, ${TGTS})

xxxx:
	@echo THIS MAKEFILE NEEDS TO BE ADAPTED FURTHER AFTER MPLEX MERGE.

all: ${EXES}

matriplex-auto:
	${MAKE} -C ../Matriplex auto && touch $@

auto: matriplex-auto
.PHONY: auto

%.o: %.cxx *.h matriplex-auto
	icc ${CPPFLAGS} ${CXXFLAGS} ${VECHOST} -c -o $@ $<

%.om: %.cxx *.h matriplex-auto
	icc ${CPPFLAGS} ${CXXFLAGS} ${VECMIC} -c -o $@ $<


clean:
	rm -f ${EXES} *.o *.om 

distclean: clean

echo:
	@echo CPPFLAGS = ${CPPFLAGS}
	@echo CXXFLAGS = ${CXXFLAGS}
	@echo LDFLAGS  = ${LDFLAGS}

################################################################

### mkFit test

MKFSRCS := $(wildcard *.cc)
MKFHDRS := $(wildcard *.h)

MKFOBJS     := $(MKFSRCS:.cc=.o)  Timing.o
MKFOBJS_MIC := $(MKFSRCS:.cc=.om) Timing.om

### To run without validation
MK_HOST_CFLAGS := -DNO_ROOT
MK_HOST_LIBS   :=
### To run with validation (host only)
# MK_HOST_CFLAGS := -I${ROOTSYS}/include
# MK_HOST_LIBS   := -L${ROOTSYS}/lib -lCore -lRIO -lTree

mkFit: auto ${MKFOBJS}
	icc ${CXXFLAGS} ${VECHOST} ${LDFLAGS} ${MK_HOST_LIBS} -o $@ ${MKFOBJS}

mkFit-mic: auto ${MKFOBJS_MIC}
	icc ${CXXFLAGS} ${VECMIC}  ${LDFLAGS} -o $@ ${MKFOBJS_MIC}
	scp $@ mic0:

%.o: %.cc *.h ../Matriplex/*
	icc ${CPPFLAGS} ${CXXFLAGS} ${VECHOST} -IMatriplex ${MK_HOST_CFLAGS} -c -o $@ $<

%.om: %.cc *.h ../Matriplex/*
	icc ${CPPFLAGS} ${CXXFLAGS} ${VECMIC} -DNO_ROOT -IMatriplex -c -o $@ $<
