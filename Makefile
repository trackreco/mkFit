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

distclean: clean
	rm -f ${AUTO_TGTS}


main: ${AUTO_TGTS} ${OBJS}
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${LDFLAGS} -o $@ ${OBJS}

${OBJS}: %.o: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<


ifeq ($(CXX),icc)

OBJS_MIC := $(OBJS:.o=.om)

main-mic: ${AUTO_TGTS} ${OBJS_MIC}
	${CXX} ${CXXFLAGS} ${VEC_MIC} ${LDFLAGS} -o $@ ${OBJS_MIC}
	scp $@ mic0:

${OBJS_MIC}: %.om: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_MIC} -DNO_ROOT -c -o $@ $<

endif
