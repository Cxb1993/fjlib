CC_LINKER=mpiCC
#GPROF=-g
GPROF=-DNDEBUG -O2
CFLAGS=
FFLAGS=
CPPFLAGS=-w -I../../include ${GPROF}
PPPFLAGS=

include ${PETSC_DIR}/bmake/common/base

SOURCE=bubble.cpp fjapp_ksp.cpp fjapp_kspm.cpp fjapp_bubbleNS.cpp fjapp_bubbleNSG.cpp fjapp_bubbleNSGD.cpp fjapp_SurfTsnEQ.cpp fjlib_cubicspline.cpp fjlib_curve.cpp fjlib_interpolation.cpp fjlib_memini.cpp fjlib_polynomial.cpp fjlib_polyroots.cpp fjlib_smoothing.cpp fjlib_solver.cpp fjlib_root.cpp

#EXTRA_DEP=fjapp_bubble_driver.cpp fjapp_vof.cpp fjapp_vof_staggered.cpp

vpath %.cpp ../../src

PROJECT=bubble

OBJECTS=${SOURCE:.cpp=.o}

${PROJECT}: ${OBJECTS} chkopts

	-${CLINKER} -o $@ ${OBJECTS} ${PETSC_KSP_LIB} ${GPROF}

all: ${PROJECT}

