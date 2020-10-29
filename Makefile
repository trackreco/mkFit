include Makefile.config

LIB_CORE     := lib/libMicCore.so
LIB_CORE_MIC := lib/libMicCore-mic.so


TGTS := ${LIB_CORE} main

ifdef KNC_BUILD
  TGTS += ${LIB_CORE_MIC} main-mic
endif

.PHONY: all clean distclean

all: ${TGTS}
	cd Geoms && ${MAKE}
	cd mkFit && ${MAKE}

SRCS := $(wildcard *.cc)
OBJS := $(SRCS:.cc=.o)
DEPS := $(SRCS:.cc=.d)

CORE_OBJS := $(filter-out main.o, ${OBJS})

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
	-rm -f ${TGTS} *.d *.o *.om *.so
	-rm -rf main.dSYM
	-rm -rf plotting/*.so plotting/*.d plotting/*.pcm

clean: clean-local
	cd mkFit && ${MAKE} clean

distclean: clean-local
	-rm -f ${AUTO_TGTS}
	-rm -f *.optrpt
	-rm -rf lib
	cd Geoms     && ${MAKE} distclean
	cd Matriplex && ${MAKE} distclean
	cd mkFit     && ${MAKE} distclean

${LIB_CORE}: ${CORE_OBJS}
	@mkdir -p $(@D)
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${CORE_OBJS} -shared -o $@ ${LDFLAGS_HOST} ${LDFLAGS}

main: ${AUTO_TGTS} ${LIB_CORE} main.o
	${CXX} ${CXXFLAGS} ${VEC_HOST} -o $@ main.o ${LDFLAGS_HOST} ${LDFLAGS} -Llib -lMicCore -Wl,-rpath,lib

${OBJS}: %.o: %.cc %.d
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<


ifdef KNC_BUILD

OBJS_MIC      := $(OBJS:.o=.om)
CORE_OBJS_MIC := $(CORE_OBJS:.o=.om)

${LIB_CORE_MIC}: ${CORE_OBJS_MIC}
	${CXX} ${CXXFLAGS} ${VEC_MIC} ${LDFLAGS_NO_ROOT} ${CORE_OBJS_MIC} -shared -o $@ ${LDFLAGS_MIC}

main-mic: ${AUTO_TGTS} ${LIB_CORE_MIC} main.om 
	${CXX} ${CXXFLAGS} ${VEC_MIC} ${LDFLAGS_NO_ROOT} -o $@ main.om ${LDFLAGS_MIC} -Llib -lMicCore-mic -Wl,-rpath=lib

${OBJS_MIC}: %.om: %.cc
	${CXX} ${CPPFLAGS_NO_ROOT} ${CXXFLAGS} ${VEC_MIC} -c -o $@ $<

endif


echo:
	-echo CXX = ${CXX}

echo_cc_defs:
	${CXX} -dM -E -mavx2 - < /dev/null
