include Makefile.config

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
