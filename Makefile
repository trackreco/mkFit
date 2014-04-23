mult66: mult66.cc Matriplex.h Makefile
	icc -O3 -mavx -I. -o mult66 mult66.cc
	icc -O3 -mmic -I. -o mult66-mic mult66.cc
	scp mult66-mic root@mic0:

mult66_test:
	./mult66
	ssh root@mic0 ./mult66-mic
