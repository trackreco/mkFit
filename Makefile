# Requires some latest gcc, e.g.:
# . /opt/rh/devtoolset-2/enable

MPLEXDEFS := -I. -DMDIM=6
LDFLAGS :=  $(shell root-config --libs)
#CXX=icc

ifeq ($(CXX),c++)
	CXXFLAGS := -std=c++11 -O3 -openmp -Wall -Wno-unknown-pragmas -I. $(shell root-config --cflags)
	MPLEXOPTS := -std=c++11 -O3 -openmp
else
	CXX := icc
	CXXFLAGS := -std=gnu++0x -O3 -openmp -I. $(shell root-config --cflags)
	MPLEXOPTS := -std=gnu++0x -O3 -openmp -vec-report=1 # -vec-threshold=0
endif

MOBJ = main.o Matrix.o KalmanUtils.o Propagation.o Simulation.o buildtest.o fittest.o

.PHONY: all clean 

all: main

clean:
	-rm -f main *.o

main: $(MOBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

main.o: fittest.h buildtest.h
Matrix.o: Matrix.h
KalmanUtils.o: KalmanUtils.h
Propagation.o: Propagation.h
Simulation.o: Simulation.h
buildtest.o: buildtest.h KalmanUtils.h Simulation.h
fittest.o: fittest.h KalmanUtils.h Simulation.h

Hit.h: Matrix.h
KalmanUtils.h: Track.h
Propagation.h: Track.h
Simulation.h: Propagation.h
Track.h: Hit.h Matrix.h
