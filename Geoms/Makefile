include ../Makefile.config

CPPEXTRA := ${USER_CPPFLAGS} ${DEFS} -I..
LDEXTRA  := ${USER_LDFLAGS}

CPPFLAGS := ${CPPEXTRA} ${CPPFLAGS}
CXXFLAGS += -fPIC ${USER_CXXFLAGS}
LDFLAGS  += ${LDEXTRA}

.PHONY: all clean distclean echo

SRCS := $(wildcard *.cc)
DEPS := $(SRCS:.cc=.d)

TGTS := $(SRCS:.cc=)

ifeq (${CXX},icc)
  TGTS += $(addsuffix -mic, ${TGTS})
endif

TGTS := $(addsuffix .so, ${TGTS})

all: ${TGTS}

%.o: %.cc %.d
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<

%.so: %.o
	${CXX} -shared -o $@ $<

ifeq ($(filter clean distclean, ${MAKECMDGOALS}),)
include ${DEPS}
endif

clean:
	rm -f *.so *.o *.d
