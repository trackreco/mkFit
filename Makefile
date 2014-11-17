include Makefile.config

TGTS := main

EXES := ${TGTS}

ifeq (${CXX},icc)
  EXES   += $(addsuffix -mic, ${TGTS})
endif

all: ${EXES}

AUTO_TGTS :=

ifdef USE_MATRIPLEX

auto-matriplex:
	${MAKE} -C Matriplex auto && touch $@

AUTO_TGTS += auto-matriplex

endif

SRCS := $(wildcard *.cc)
DEPS := $(SRCS:.cc=.d)
OBJS := $(SRCS:.cc=.o)

-include ${DEPS}

.PHONY: all clean 

clean:
	-rm -f ${EXES} *.d *.o *.om
	-rm -rf USolids-{host,mic}

distclean: clean
	rm -f ${AUTO_TGTS}


main: ${AUTO_TGTS} ${OBJS} ${LIBUSOLIDS}
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${LDFLAGS} -o $@ ${OBJS} ${LIBUSOLIDS}

${OBJS}: %.o: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<

${LIBUSOLIDS} : USolids/CMakeLists.txt
	-mkdir USolids-host
	cd USolids-host && cmake ${CMAKEFLAGS} ../USolids && make

ifeq ($(CXX),icc)

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
