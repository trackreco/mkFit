# Requires some latest gcc, e.g.:
# . /opt/rh/devtoolset-2/enable

MPLEXDEFS := -DMDIM=3 -I.
MPLEXOPTS := -std=gnu++0x -O3 -openmp #-vec-report=2 -vec-threshold=0

%: %.cxx mplex-common.h Matriplex.h MatriplexVector.h MatriplexNT.h Makefile
	icc ${MPLEXDEFS} ${MPLEXOPTS} -mavx -o $@ $<

%-mic: %.cxx mplex-common.h Matriplex.h MatriplexVector.h MatriplexNT.h Makefile
	icc ${MPLEXDEFS} ${MPLEXOPTS} -mmic -o $@ $<
	scp $@ root@mic0:


%-test: % %-mic
	./$*
	ssh root@mic0 ./$*-mic
