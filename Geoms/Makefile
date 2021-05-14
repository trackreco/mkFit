include ../Makefile.config

CPPEXTRA := ${USER_CPPFLAGS} ${DEFS} -I..
LDEXTRA  := ${USER_LDFLAGS}

CPPFLAGS := ${CPPEXTRA} ${CPPFLAGS}
CXXFLAGS += -fPIC ${USER_CXXFLAGS}
LDFLAGS  += ${LDEXTRA}

CPPFLAGS_NO_ROOT := ${CPPEXTRA} ${CPPFLAGS_NO_ROOT}

.PHONY: all clean distclean echo

SRCS := $(wildcard *.cc)
DEPS := $(SRCS:.cc=.d)

TGTS := $(SRCS:.cc=)

TGTS := $(addsuffix .so, ${TGTS})

all: ${TGTS}

%.o: %.cc %.d
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<

%.so: %.o
	${CXX} -shared -L../lib -lMicCore -o $@ $<

ifeq ($(filter clean distclean, ${MAKECMDGOALS}),)
include ${DEPS}
endif

clean:
	rm -f *.so *.o *.om *.d *.optrpt

distclean: clean
