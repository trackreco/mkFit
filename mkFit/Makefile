##### Random stuff #####
#
# To disable vectorization set USER_CXXFLAGS := -no-simd -no-vec
# Setting only one of the above has little effect.
# Note, this also screws-up prefetching so it's a lousy deal.
#
# -opt-prefetch-distance=64,8

include ../Makefile.config

CPPEXTRA := -I.. ${USER_CPPFLAGS} ${DEFS}
LDEXTRA  := ${USER_LDFLAGS}

CPPFLAGS += ${CPPEXTRA}
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

auto-matriplex:
	${MAKE} -C ../Matriplex auto && touch $@

AUTO_HEADERS := upParam_MultKalmanGain.ah \
		upParam_simil_x_propErr.ah \
		upParam_propErrT_x_simil_propErr.ah \
		upParam_kalmanGain_x_propErr.ah \
		MultHelixProp.ah \
		MultHelixPropTransp.ah

${AUTO_HEADERS}: auto-genmplex

.INTERMEDIATE: auto-genmplex
auto-genmplex: GenMPlexOps.pl
	./GenMPlexOps.pl

AUTO_TGTS := auto-matriplex ${AUTO_HEADERS}

default: ${AUTO_TGTS} ${EXES}

clean:
	rm -f ${EXES} *.d *.o *.om 

distclean: clean
	rm -f *.optrpt
	rm -f ${AUTO_TGTS}

echo:
	@echo "CXX      = ${CXX}"
	@echo "CPPFLAGS = ${CPPFLAGS}"
	@echo "CXXFLAGS = ${CXXFLAGS}"
	@echo "LDFLAGS  = ${LDFLAGS}"
	@echo "EXES     = ${EXES}"


################################################################

# Should be a lib, really
ABOVE_OBJS := $(patsubst %, ../%.o, Config Matrix Event Hit Track Propagation KalmanUtils Simulation Geometry SimpleGeom fittest buildtest ConformalUtils seedtest BinInfoUtils)

${ABOVE_OBJS}:
	${MAKE} -C ..

MKFSRCS := $(wildcard *.cc)
MKFHDRS := $(wildcard *.h)

MKFOBJS := $(MKFSRCS:.cc=.o)

ALLOBJS := ${MKFOBJS} ${ABOVE_OBJS}

ifeq ($(filter clean distclean, ${MAKECMDGOALS}),)
include $(MKFSRCS:.cc=.d)
endif

mkFit: ${ALLOBJS}
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${LDFLAGS} ${ALLOBJS} -o $@ ${LDFLAGS_HOST}

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
