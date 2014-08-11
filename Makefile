# Requires some latest gcc, e.g.:
# . /opt/rh/devtoolset-2/enable

MPLEXDEFS := -I. -DMDIM=6
LDFLAGS :=  $(shell . /Volumes/User/Applications/root/bin/thisroot.sh && root-config --libs)
#CXX=icc

ifeq ($(CXX),c++)
	CXXFLAGS := -std=c++11 -O3 -openmp -Wall -Wno-unknown-pragmas -I. $(shell . /Volumes/User/Applications/root/bin/thisroot.sh && root-config --cflags)
	MPLEXOPTS := -std=c++11 -O3 -openmp
else
	CXX := icc
	CXXFLAGS := -std=gnu++0x -O3 -openmp -I. $(shell . /Volumes/User/Applications/root/bin/thisroot.sh && root-config --cflags)
	MPLEXOPTS := -std=gnu++0x -O3 -openmp -vec-report=1 # -vec-threshold=0
endif

MOBJ = main.o Matrix.o KalmanUtils.o Propagation.o Simulation.o buildtest.o fittest.o Geometry.o

MEXE-MIC = mplex-mic mplex-vec-mic mplex-nt-mic mplexsym-mic mplexsym-nt-mic
MEXE-AVX = mplex mplex-vec mplex-nt mplexsym mplexsym-nt

all: main #all-avx all-mic

all-mic:	$(MEXE-MIC)
all-avx:	$(MEXE-AVX)

test:	mplex-test mplex-vec-test mplex-nt-test mplexsym-test mplexsym-nt-test

%: %.cxx mplex-common.h Matriplex.h MatriplexSym.h MatriplexVector.h MatriplexNT.h MatriplexSymNT.h Makefile
	$(CXX) ${MPLEXDEFS} ${MPLEXOPTS} -mavx -o $@ $< mplex-common.cxx

%-mic: %.cxx mplex-common.h Matriplex.h MatriplexSym.h MatriplexVector.h MatriplexNT.h MatriplexSymNT.h Makefile
	$(CXX) ${MPLEXDEFS} ${MPLEXOPTS} -mmic -o $@ $< mplex-common.cxx

%-test: % %-mic
	./$*
	scp $*-mic root@mic0:
	ssh root@mic0 ./$*-mic

main: $(MOBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) -LUSolids -lusolids

buildtest.o: buildtest.cc Propagation.h buildtest.h KalmanUtils.h Simulation.h Geometry.h
fittest.o: fittest.cc Propagation.h fittest.h KalmanUtils.h Simulation.h Geometry.h
Geometry.o: Geometry.h
KalmanUtils.o: KalmanUtils.cc KalmanUtils.h Track.h
main.o: main.cc fittest.h buildtest.h
Matrix.o: Matrix.cc Matrix.h
KalmanUtils.o: KalmanUtils.cc KalmanUtils.h
Propagation.o: Propagation.cc Propagation.h Track.h
Simulation.o: Simulation.cc Propagation.h Simulation.h Geometry.h

Hit.h: Matrix.h
KalmanUtils.h: Track.h
Matriplex.h: MatriplexCommon.h
MatriplexNT.h: MatriplexCommon.h
MatriplexSym.h: MatriplexCommon.h Matriplex.h
MatriplexSymNT.h: MatriplexCommon.h MatriplexNT.h
MatriplexVector.h: Matriplex.h
Matrix.h: MatriplexSymNT.h
Propagation.h: Track.h
Simulation.h: Propagation.h
Track.h: Hit.h Matrix.h

.PHONY: all all-mic all-avx test clean 

clean:
	-rm $(MEXE-AVX) $(MEXE-MIC)
