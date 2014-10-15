include Makefile.config

TGTS := main

EXES := ${TGTS}

ifeq (${CXX},icc)
  EXES   += $(addsuffix -mic, ${TGTS})
endif


all: ${EXES}


SRCS := $(wildcard *.cc)
DEPS := $(SRCS:.cc=.d)
OBJS := $(SRCS:.cc=.o)

-include ${DEPS}

.PHONY: all clean 

clean:
	-rm -f ${EXES} *.d *.o *.om

distclean: clean


main: ${OBJS}
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${LDFLAGS} -o $@ $^

${OBJS}: %.o: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<


ifeq ($(CXX),icc)

OBJS_MIC := $(OBJS:.o=.om)

main-mic: ${OBJS_MIC}
	${CXX} ${CXXFLAGS} ${VEC_MIC} ${LDFLAGS} -o $@ $^
	scp $@ mic0:

${OBJS_MIC}: %.om: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_MIC} -DNO_ROOT -c -o $@ $<

endif
