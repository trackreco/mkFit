
.PHONY: all clean distclean echo

ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else
all: default

CPPEXTRA := -I.. ${USER_CPPFLAGS} ${DEFS}
LDEXTRA  := ${USER_LDFLAGS}

CPPFLAGS := $(shell root-config --cflags) ${CPPEXTRA} ${CPPFLAGS}
CXXFLAGS += ${USER_CXXFLAGS}
LDFLAGS  += $(shell root-config --ldflags --cflags --libs) ${LDEXTRA}


TGTS     := writeMemoryFile

EXES     := ${TGTS}

default: ${EXES}

clean:
	rm -f ${EXES} *.d *.o *.om
	rm -rf writeMemoryFile.dSYM

distclean: clean
	rm -f *.optrpt 
	rm -f ${EXES} *.ah

echo:
	@echo "CXX      = ${CXX}"
	@echo "CPPFLAGS = ${CPPFLAGS}"
	@echo "CXXFLAGS = ${CXXFLAGS}"
	@echo "LDFLAGS  = ${LDFLAGS}"
	@echo "EXES     = ${EXES}"


################################################################

# Should be a lib, really
ABOVE_OBJS := $(patsubst %, ../%.o, Config Matrix Hit Track )

${ABOVE_OBJS}:
	${MAKE} -C ..

SRCS := $(wildcard *.cc)
HDRS := $(wildcard *.h)

OBJS := $(SRCS:.cc=.o)

ALLOBJS := ${OBJS} ${ABOVE_OBJS}


writeMemoryFile: ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${ALLOBJS} -o $@


${OBJS}: %.o: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c -o $@ $<
endif

