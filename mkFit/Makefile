##### Random stuff #####
#
# To disable vectorization set USER_CXXFLAGS := -no-simd -no-vec
# Setting only one of the above has little effect.
# Note, this also screws-up prefetching so it's a lousy deal.
#
# -opt-prefetch-distance=64,8

include ../Makefile.config


CPPFLAGS += -I.. ${USER_CPPFLAGS} ${DEFS}
CXXFLAGS += ${USER_CXXFLAGS}

LDFLAGS  += ${USER_LDFLAGS}

TGTS     := mkFit

EXES     := ${TGTS}

ifeq (${CXX},icc)
  EXES   += $(addsuffix -mic, ${TGTS})
endif


all: ${EXES}


auto-matriplex:
	${MAKE} -C ../Matriplex auto && touch $@

auto-genmplex: GenMPlexOps.pl
	./GenMPlexOps.pl && touch $@

AUTO_TGTS := auto-matriplex auto-genmplex


clean:
	rm -f ${EXES} *.d *.o *.om 

distclean: clean
	rm -f ${AUTO_TGTS} *.ah

echo:
	@echo "CXX      = ${CXX}"
	@echo "CPPFLAGS = ${CPPFLAGS}"
	@echo "CXXFLAGS = ${CXXFLAGS}"
	@echo "LDFLAGS  = ${LDFLAGS}"
	@echo "EXES     = ${EXES}"


################################################################

-include *.d ../*.d

# Should be a lib, really
ABOVE_OBJS := $(patsubst %, ../%.o, Matrix Hit Track Propagation KalmanUtils Simulation)

${ABOVE_OBJS}:
	${MAKE} -C ..

MKFSRCS := $(wildcard *.cc)
MKFHDRS := $(wildcard *.h)

MKFOBJS := $(MKFSRCS:.cc=.o)

ALLOBJS := ${MKFOBJS} ${ABOVE_OBJS}

mkFit: ${AUTO_TGTS} ${ALLOBJS}
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${LDFLAGS} -o $@ ${ALLOBJS}

${MKFOBJS}: %.o: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<


### MIC build, icc only

ifeq ($(CXX),icc)

MKFOBJS_MIC := $(MKFOBJS:.o=.om)

ALLOBJS_MIC := $(ALLOBJS:.o=.om) 

mkFit-mic: ${AUTO_TGTS} ${ALLOBJS_MIC}
	${CXX} ${CXXFLAGS} ${VEC_MIC} ${LDFLAGS} -o $@ ${ALLOBJS_MIC}
	scp $@ mic0:

${MKFOBJS_MIC}: %.om: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_MIC} -DNO_ROOT -IMatriplex -c -o $@ $<

endif
