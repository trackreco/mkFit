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

MOBJ = main.o Matrix.o KalmanUtils.o Propagation.o Simulation.o buildtest.o fittest.o Geometry.o ConformalUtils.o

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
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) -LUSolids -lusolids

.PHONY: all all-mic all-avx test clean depend

clean:
	-rm $(MEXE-AVX) $(MEXE-MIC) *.o
depend:
	makedepend -- $(CFLAGS) --  *.cc *.h

# DO NOT DELETE

ConformalUtils.o: ConformalUtils.h Track.h Hit.h Matrix.h Math/SMatrix.h
ConformalUtils.o: Math/MConfig.h Math/Expression.h
ConformalUtils.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
ConformalUtils.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
ConformalUtils.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
ConformalUtils.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
ConformalUtils.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
ConformalUtils.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h
ConformalUtils.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
ConformalUtils.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h
Geometry.o: Geometry.h USolids/include/VUSolid.hh USolids/include/UTypes.hh
Geometry.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
Hit.o: Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h Math/Expression.h
Hit.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h Math/SMatrix.icc
Hit.o: Math/Dfact.h Math/Dinv.h Math/CholeskyDecomp.h
Hit.o: Math/CramerInversion.icc Math/CramerInversionSym.icc
Hit.o: Math/MatrixInversion.icc Math/SVector.h Math/SVector.icc
Hit.o: Math/UnaryOperators.h Math/BinaryOperators.h Math/BinaryOpPolicy.h
Hit.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
Hit.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h
KalmanUtils.o: KalmanUtils.h Track.h Hit.h Matrix.h Math/SMatrix.h
KalmanUtils.o: Math/MConfig.h Math/Expression.h
KalmanUtils.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
KalmanUtils.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
KalmanUtils.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
KalmanUtils.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
KalmanUtils.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
KalmanUtils.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h Math/Functions.h
KalmanUtils.o: Math/HelperOps.h Math/MatrixFunctions.h MatriplexSymNT.h
KalmanUtils.o: MatriplexCommon.h MatriplexNT.h KalmanOpsNT.h
Matrix.o: Matrix.h Math/SMatrix.h Math/MConfig.h Math/Expression.h
Matrix.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
Matrix.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h Math/CholeskyDecomp.h
Matrix.o: Math/CramerInversion.icc Math/CramerInversionSym.icc
Matrix.o: Math/MatrixInversion.icc Math/SVector.h Math/SVector.icc
Matrix.o: Math/UnaryOperators.h Math/BinaryOperators.h Math/BinaryOpPolicy.h
Matrix.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
Matrix.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h MatriplexSymNT.icc
Propagation.o: Propagation.h Track.h Hit.h Matrix.h Math/SMatrix.h
Propagation.o: Math/MConfig.h Math/Expression.h
Propagation.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
Propagation.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
Propagation.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
Propagation.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
Propagation.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
Propagation.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h Math/Functions.h
Propagation.o: Math/HelperOps.h Math/MatrixFunctions.h MatriplexSymNT.h
Propagation.o: MatriplexCommon.h MatriplexNT.h Geometry.h
Propagation.o: USolids/include/VUSolid.hh USolids/include/UTypes.hh
Propagation.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
Simulation.o: Simulation.h Track.h Hit.h Matrix.h Math/SMatrix.h
Simulation.o: Math/MConfig.h Math/Expression.h
Simulation.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
Simulation.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h Math/CholeskyDecomp.h
Simulation.o: Math/CramerInversion.icc Math/CramerInversionSym.icc
Simulation.o: Math/MatrixInversion.icc Math/SVector.h Math/SVector.icc
Simulation.o: Math/UnaryOperators.h Math/BinaryOperators.h
Simulation.o: Math/BinaryOpPolicy.h Math/Functions.h Math/HelperOps.h
Simulation.o: Math/MatrixFunctions.h MatriplexSymNT.h MatriplexCommon.h
Simulation.o: MatriplexNT.h Propagation.h Geometry.h
Simulation.o: USolids/include/VUSolid.hh USolids/include/UTypes.hh
Simulation.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
Track.o: Track.h Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h
Track.o: Math/Expression.h Math/MatrixRepresentationsStatic.h
Track.o: Math/StaticCheck.h Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
Track.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
Track.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc Math/SVector.h
Track.o: Math/SVector.icc Math/UnaryOperators.h Math/BinaryOperators.h
Track.o: Math/BinaryOpPolicy.h Math/Functions.h Math/HelperOps.h
Track.o: Math/MatrixFunctions.h MatriplexSymNT.h MatriplexCommon.h
Track.o: MatriplexNT.h
buildtest.o: buildtest.h Track.h Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h
buildtest.o: Math/Expression.h Math/MatrixRepresentationsStatic.h
buildtest.o: Math/StaticCheck.h Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
buildtest.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
buildtest.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
buildtest.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
buildtest.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h Math/Functions.h
buildtest.o: Math/HelperOps.h Math/MatrixFunctions.h MatriplexSymNT.h
buildtest.o: MatriplexCommon.h MatriplexNT.h Geometry.h
buildtest.o: USolids/include/VUSolid.hh USolids/include/UTypes.hh
buildtest.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
buildtest.o: KalmanUtils.h Propagation.h Simulation.h
fittest.o: fittest.h Geometry.h USolids/include/VUSolid.hh
fittest.o: USolids/include/UTypes.hh USolids/include/UVector3.hh
fittest.o: USolids/include/UUtils.hh Track.h Hit.h Matrix.h Math/SMatrix.h
fittest.o: Math/MConfig.h Math/Expression.h
fittest.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
fittest.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h Math/CholeskyDecomp.h
fittest.o: Math/CramerInversion.icc Math/CramerInversionSym.icc
fittest.o: Math/MatrixInversion.icc Math/SVector.h Math/SVector.icc
fittest.o: Math/UnaryOperators.h Math/BinaryOperators.h Math/BinaryOpPolicy.h
fittest.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
fittest.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h KalmanUtils.h
fittest.o: Propagation.h Simulation.h ConformalUtils.h
main.o: Matrix.h Math/SMatrix.h Math/MConfig.h Math/Expression.h
main.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
main.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h Math/CholeskyDecomp.h
main.o: Math/CramerInversion.icc Math/CramerInversionSym.icc
main.o: Math/MatrixInversion.icc Math/SVector.h Math/SVector.icc
main.o: Math/UnaryOperators.h Math/BinaryOperators.h Math/BinaryOpPolicy.h
main.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
main.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h fittest.h Geometry.h
main.o: USolids/include/VUSolid.hh USolids/include/UTypes.hh
main.o: USolids/include/UVector3.hh USolids/include/UUtils.hh buildtest.h
main.o: Track.h Hit.h USolids/include/UTubs.hh USolids/include/VUSolid.hh
main.o: USolids/include/UTubs.icc USolids/include/UPolyHedra.hh
main.o: USolids/include/UVCSGfaceted.hh USolids/include/UVoxelizer.hh
main.o: USolids/include/UBits.hh USolids/include/UBox.hh
main.o: USolids/include/VUFacet.hh USolids/include/UTransform3D.hh
main.o: USolids/include/UReduciblePolygon.hh
main.o: USolids/include/UPolyhedraSide.hh USolids/include/UVCSGface.hh
main.o: USolids/include/UPolyhedra.icc USolids/include/USphere.hh
ConformalUtils.o: Track.h Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h
ConformalUtils.o: Math/Expression.h Math/MatrixRepresentationsStatic.h
ConformalUtils.o: Math/StaticCheck.h Math/SMatrix.icc Math/Dfact.h
ConformalUtils.o: Math/Dinv.h Math/CholeskyDecomp.h Math/CramerInversion.icc
ConformalUtils.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
ConformalUtils.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
ConformalUtils.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h
ConformalUtils.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
ConformalUtils.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h
Geometry.o: USolids/include/VUSolid.hh USolids/include/UTypes.hh
Geometry.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
Hit.o: Matrix.h Math/SMatrix.h Math/MConfig.h Math/Expression.h
Hit.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h Math/SMatrix.icc
Hit.o: Math/Dfact.h Math/Dinv.h Math/CholeskyDecomp.h
Hit.o: Math/CramerInversion.icc Math/CramerInversionSym.icc
Hit.o: Math/MatrixInversion.icc Math/SVector.h Math/SVector.icc
Hit.o: Math/UnaryOperators.h Math/BinaryOperators.h Math/BinaryOpPolicy.h
Hit.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
Hit.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h
KalmanUtils.o: Track.h Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h
KalmanUtils.o: Math/Expression.h Math/MatrixRepresentationsStatic.h
KalmanUtils.o: Math/StaticCheck.h Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
KalmanUtils.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
KalmanUtils.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
KalmanUtils.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
KalmanUtils.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h Math/Functions.h
KalmanUtils.o: Math/HelperOps.h Math/MatrixFunctions.h MatriplexSymNT.h
KalmanUtils.o: MatriplexCommon.h MatriplexNT.h
Matriplex.o: MatriplexCommon.h
MatriplexNT.o: MatriplexCommon.h
MatriplexSym.o: MatriplexCommon.h Matriplex.h
MatriplexSymNT.o: MatriplexCommon.h MatriplexNT.h
MatriplexVector.o: Matriplex.h MatriplexCommon.h
Matrix.o: Math/SMatrix.h Math/MConfig.h Math/Expression.h
Matrix.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
Matrix.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h Math/CholeskyDecomp.h
Matrix.o: Math/CramerInversion.icc Math/CramerInversionSym.icc
Matrix.o: Math/MatrixInversion.icc Math/SVector.h Math/SVector.icc
Matrix.o: Math/UnaryOperators.h Math/BinaryOperators.h Math/BinaryOpPolicy.h
Matrix.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
Matrix.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h
Propagation.o: Track.h Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h
Propagation.o: Math/Expression.h Math/MatrixRepresentationsStatic.h
Propagation.o: Math/StaticCheck.h Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
Propagation.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
Propagation.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
Propagation.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
Propagation.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h Math/Functions.h
Propagation.o: Math/HelperOps.h Math/MatrixFunctions.h MatriplexSymNT.h
Propagation.o: MatriplexCommon.h MatriplexNT.h Geometry.h
Propagation.o: USolids/include/VUSolid.hh USolids/include/UTypes.hh
Propagation.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
Simulation.o: Track.h Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h
Simulation.o: Math/Expression.h Math/MatrixRepresentationsStatic.h
Simulation.o: Math/StaticCheck.h Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
Simulation.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
Simulation.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
Simulation.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
Simulation.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h Math/Functions.h
Simulation.o: Math/HelperOps.h Math/MatrixFunctions.h MatriplexSymNT.h
Simulation.o: MatriplexCommon.h MatriplexNT.h Propagation.h Geometry.h
Simulation.o: USolids/include/VUSolid.hh USolids/include/UTypes.hh
Simulation.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
Track.o: Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h Math/Expression.h
Track.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
Track.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h Math/CholeskyDecomp.h
Track.o: Math/CramerInversion.icc Math/CramerInversionSym.icc
Track.o: Math/MatrixInversion.icc Math/SVector.h Math/SVector.icc
Track.o: Math/UnaryOperators.h Math/BinaryOperators.h Math/BinaryOpPolicy.h
Track.o: Math/Functions.h Math/HelperOps.h Math/MatrixFunctions.h
Track.o: MatriplexSymNT.h MatriplexCommon.h MatriplexNT.h
buildtest.o: Track.h Hit.h Matrix.h Math/SMatrix.h Math/MConfig.h
buildtest.o: Math/Expression.h Math/MatrixRepresentationsStatic.h
buildtest.o: Math/StaticCheck.h Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
buildtest.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
buildtest.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
buildtest.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
buildtest.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h Math/Functions.h
buildtest.o: Math/HelperOps.h Math/MatrixFunctions.h MatriplexSymNT.h
buildtest.o: MatriplexCommon.h MatriplexNT.h Geometry.h
buildtest.o: USolids/include/VUSolid.hh USolids/include/UTypes.hh
buildtest.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
fittest.o: Geometry.h USolids/include/VUSolid.hh USolids/include/UTypes.hh
fittest.o: USolids/include/UVector3.hh USolids/include/UUtils.hh
mplex-common.o: Math/SMatrix.h Math/MConfig.h Math/Expression.h
mplex-common.o: Math/MatrixRepresentationsStatic.h Math/StaticCheck.h
mplex-common.o: Math/SMatrix.icc Math/Dfact.h Math/Dinv.h
mplex-common.o: Math/CholeskyDecomp.h Math/CramerInversion.icc
mplex-common.o: Math/CramerInversionSym.icc Math/MatrixInversion.icc
mplex-common.o: Math/SVector.h Math/SVector.icc Math/UnaryOperators.h
mplex-common.o: Math/BinaryOperators.h Math/BinaryOpPolicy.h Math/Functions.h
mplex-common.o: Math/HelperOps.h Math/MatrixFunctions.h
