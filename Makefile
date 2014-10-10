include Makefile.config

CPPFLAGS := -I. -MMD
CXXFLAGS := -O3 -std=c++11
LDFLAGS  :=

ifeq (${WITH_ROOT},yes)
  CPPFLAGS += $(shell root-config --cflags)
  CXXFLAGS += 
  LDFLAGS  += $(shell root-config --libs)
else
  CPPFLAGS += -DNO_ROOT
endif

ifeq ($(CXX),icc)
  CXXFLAGS += -openmp
  LDFLAGS  += -openmp
else
  CXXFLAGS += -fopenmp -Wall -Wno-unknown-pragmas
  LDFLAGS  += -fopenmp
endif

all: main

SRCS := $(wildcard *.cc)
DEPS := $(SRCS:.cc=.d)
OBJS := $(SRCS:.cc=.o)

-include ${DEPS}

.PHONY: all clean 

clean:
	-rm -f main *.o *.d

main: ${OBJS}
	${CXX} -o $@ $^ ${LDFLAGS}
