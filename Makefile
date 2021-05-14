include Makefile.config

LIB_CORE     := lib/libMicCore.so

TGTS := ${LIB_CORE}

.PHONY: all clean distclean

all: ${TGTS}
	cd Geoms && ${MAKE}
	cd mkFit && ${MAKE}

SRCS := $(wildcard *.cc)
OBJS := $(SRCS:.cc=.o)
DEPS := $(SRCS:.cc=.d)

CORE_OBJS := $(filter-out main.o, ${OBJS})

AUTO_TGTS :=

# Begin Matriplex

auto-matriplex:
	${MAKE} -C Matriplex auto && touch $@

AUTO_TGTS += auto-matriplex

${DEPS}: auto-matriplex

# End Matriplex

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

echo:
	@echo "CXX=${CXX}"

echo_cc_defs:
	${CXX} -dM -E -mavx2 - < /dev/null

echo-srcs:
	@echo ${SRCS}

echo-flags:
	@echo "CPPFLAGS=${CPPFLAGS}"

echo-tbb:
	@echo "TBB_GCC=${TBB_GCC}, TBB_PREFIX=${TBB_PREFIX}, TBB_ROOT=${TBB_ROOT}"
