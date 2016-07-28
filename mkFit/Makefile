##### Random stuff #####
#
# To disable vectorization set USER_CXXFLAGS := -no-simd -no-vec
# Setting only one of the above has little effect.
# Note, this also screws-up prefetching so it's a lousy deal.
#
# -opt-prefetch-distance=64,8

include ../Makefile.config

CPPEXTRA := -I.. ${USER_CPPFLAGS} ${DEFS}
LDEXTRA  := -ltbb ${USER_LDFLAGS}

CPPFLAGS := ${CPPEXTRA} ${CPPFLAGS}
CXXFLAGS += ${USER_CXXFLAGS}
LDFLAGS  += ${LDEXTRA}

CPPFLAGS_NO_ROOT += ${CPPEXTRA}
LDFLAGS_NO_ROOT  += ${LDEXTRA}

.PHONY: all clean distclean echo

all: default

TGTS     := mkFit

EXES     := ${TGTS}

ifeq (${CXX},icc)
  EXES   += $(addsuffix -mic, ${TGTS})
endif

auto-genmplex: GenMPlexOps.pl
	./GenMPlexOps.pl && touch $@

AUTO_TGTS := auto-genmplex

default: ${AUTO_TGTS} ${EXES}

clean:
	rm -f ${EXES} *.d *.o *.om Ice/*.d Ice/*.o Ice/*.om
	rm -rf mkFit.dSYM

distclean: clean
	rm -f *.optrpt Ice/*.optrpt
	rm -f ${AUTO_TGTS} *.ah

echo:
	@echo "CXX      = ${CXX}"
	@echo "CPPFLAGS = ${CPPFLAGS}"
	@echo "CXXFLAGS = ${CXXFLAGS}"
	@echo "LDFLAGS  = ${LDFLAGS}"
	@echo "EXES     = ${EXES}"
	@echo "CU_OBJ   = ${CU_OBJS}"


################################################################

# Should be a lib, really
ABOVE_OBJS := $(patsubst %, ../%.o, Config Matrix Event Hit Track Propagation KalmanUtils Simulation Geometry SimpleGeom fittest buildtest ConformalUtils seedtest BinInfoUtils TTreeValidation)

${ABOVE_OBJS}:
	${MAKE} -C ..

MKFSRCS := $(wildcard *.cc) $(wildcard Ice/*.cc)
MKFHDRS := $(wildcard *.h)

MKFOBJS := $(MKFSRCS:.cc=.o)
MKFDEPS := $(MKFSRCS:.cc=.d)

ifdef USE_CUDA
CU_SRCS := $(wildcard *.cu)
CU_OBJS := $(CU_SRCS:.cu=.o)

LDFLAGS_CU := -lcudart

# TODO: Clean the "-I.. -std=c++11"
${CU_OBJS}: %.o: %.cu
	${NV} -c -o $@ $< -I.. -std=c++11 -DUSE_MATRIPLEX
endif

ALLOBJS := ${MKFOBJS} ${ABOVE_OBJS} ${CU_OBJS}

${MKFDEPS}: auto-genmplex

ifeq ($(filter clean distclean, ${MAKECMDGOALS}),)
include ${MKFDEPS}
endif

mkFit: ${ALLOBJS}
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${LDFLAGS} ${ALLOBJS} -o $@ ${LDFLAGS_HOST} ${LDFLAGS_CU}

ifdef WITH_ROOT
fittestMPlex.o : CPPFLAGS += $(shell root-config --cflags)

TFile.h: 
	echo "Using dummy rule for TFile.h"

TTree.h: 
	echo "Using dummy rule for TTree.h"
endif

${MKFOBJS}: %.o: %.cc %.d
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<

### MIC build, icc only

ifeq ($(CXX),icc)

MKFOBJS_MIC := $(MKFOBJS:.o=.om)

ALLOBJS_MIC := $(ALLOBJS:.o=.om) 

mkFit-mic: ${ALLOBJS_MIC}
	${CXX} ${CXXFLAGS} ${VEC_MIC} ${LDFLAGS_NO_ROOT} ${ALLOBJS_MIC} -o $@ ${LDFLAGS_MIC}
	scp $@ mic0:

${MKFOBJS_MIC}: %.om: %.cc %.d
	${CXX} ${CPPFLAGS_NO_ROOT} ${CXXFLAGS} ${VEC_MIC} -c -o $@ $<

endif
