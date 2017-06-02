include ../Makefile.config

# Overide debug settings
# OPT := -O0 -g

.PHONY: all clean distclean echo

ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else
all: default

CPPEXTRA := -I.. ${USER_CPPFLAGS} ${DEFS}
LDEXTRA  := ${USER_LDFLAGS}

CPPFLAGS := $(shell root-config --cflags) ${CPPEXTRA} ${CPPFLAGS}
CXXFLAGS += -fPIC ${USER_CXXFLAGS}
LDFLAGS  += -fPIC $(shell root-config --ldflags --cflags --libs) ${LDEXTRA}


TGTS     := writeMemoryFile libDicts.so

EXES     := ${TGTS}

default: ${EXES}

clean:
	rm -f ${EXES} *.d *.o *.om *.so *Dict.* *.pcm
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

SRCS := $(wildcard *.cc)
HDRS := $(wildcard *.h)

OBJS := $(SRCS:.cc=.o)
OBJS := $(filter-out DictsDict.o, $(OBJS))

ALLOBJS := ${OBJS} ${ABOVE_OBJS}

libDicts.so: DictsDict.cc
	${CXX} $(LDFLAGS) -Wno-deprecated -shared DictsDict.cc -o $@

DictsDict.cc: DictsLinkDef.h
	@rootcint -v3 -f $@ -c -p $<

writeMemoryFile: ${OBJS} libDicts.so
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${ALLOBJS} -L. -lDicts -L.. -lMicCore -Wl,-rpath=..,-rpath=. -o $@


${OBJS}: %.o: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c -o $@ $<
endif

