include Makefile.config

TGTS := main

EXES := ${TGTS}

ifdef KNC_BUILD
  EXES   += $(addsuffix -mic, ${TGTS})
endif

.PHONY: all clean distclean

all: ${EXES}
	cd mkFit && ${MAKE}

SRCS := $(wildcard *.cc)
OBJS := $(SRCS:.cc=.o)
DEPS := $(SRCS:.cc=.d)

AUTO_TGTS :=

ifdef USE_MATRIPLEX

auto-matriplex:
	${MAKE} -C Matriplex auto && touch $@

AUTO_TGTS += auto-matriplex

${DEPS}: auto-matriplex

endif

ifeq ($(filter clean-local clean distclean, ${MAKECMDGOALS}),)
include ${DEPS}
endif

clean-local:
	-rm -f ${EXES} *.d *.o *.om *.so
	-rm -rf main.dSYM
	-rm -rf USolids-{host,mic}

clean: clean-local
	cd mkFit && ${MAKE} clean

distclean: clean-local
	-rm -f ${AUTO_TGTS}
	-rm -f *.optrpt
	cd Matriplex && ${MAKE} distclean
	cd mkFit     && ${MAKE} distclean

main: ${AUTO_TGTS} ${OBJS} ${LIBUSOLIDS}
	${CXX} ${CXXFLAGS} ${VEC_HOST} -o $@ ${OBJS} ${LIBUSOLIDS} ${LDFLAGS}

${OBJS}: %.o: %.cc %.d
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<

${LIBUSOLIDS} : USolids/CMakeLists.txt
	-mkdir USolids-host
	cd USolids-host && cmake ${CMAKEFLAGS} ../USolids && make


ifdef KNC_BUILD

OBJS_MIC := $(OBJS:.o=.om)

main-mic: ${AUTO_TGTS} ${OBJS_MIC} ${LIBUSOLIDS_MIC}
	${CXX} ${CXXFLAGS} ${VEC_MIC} ${LDFLAGS_NO_ROOT} -o $@ ${OBJS_MIC} ${LIBUSOLIDS_MIC}
	scp $@ mic0:

${LIBUSOLIDS_MIC} : USolids/CMakeLists.txt
	-mkdir USolids-mic
	cd USolids-mic && cmake ${CMAKEFLAGS_MIC} ../USolids && make

${OBJS_MIC}: %.om: %.cc
	${CXX} ${CPPFLAGS_NO_ROOT} ${CXXFLAGS} ${VEC_MIC} -c -o $@ $<

endif


echo:
	-echo CXX = ${CXX}
