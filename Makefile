# Requires some latest gcc, e.g.:
# . /opt/rh/devtoolset-2/enable

MPLEXDEFS := -I. -DMDIM=3
MPLEXOPTS := -std=gnu++0x -O3 -openmp # -vec-report=1 # -vec-threshold=0

all:	mplex mplex-vec mplex-nt mplex-mic mplex-vec-mic mplex-nt-mic mplexsym mplexsym-mic mplexsym-nt mplexsym-nt-mic

test:	mplex-test mplex-vec-test mplex-nt-test mplexsym-test mplexsym-ntst

%: %.cxx mplex-common.h Matriplex.h MatriplexSym.h MatriplexVector.h MatriplexNT.h MatriplexSymNT.h Makefile
	icc ${MPLEXDEFS} ${MPLEXOPTS} -mavx -o $@ $< mplex-common.cxx

%-mic: %.cxx mplex-common.h Matriplex.h MatriplexSym.h MatriplexVector.h MatriplexNT.h MatriplexSymNT.h Makefile
	icc ${MPLEXDEFS} ${MPLEXOPTS} -mmic -o $@ $< mplex-common.cxx
	scp $@ root@mic0:


%-test: % %-mic
	./$*
	ssh root@mic0 ./$*-mic
