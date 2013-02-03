$GPROF=-g
CC_LINKER=mpiCC
GPROF=
CFLAGS=
FFLAGS=
CPPFLAGS=-w -I../../include ${GPROF}
PPPFLAGS=

include ${PETSC_DIR}/bmake/common/base

vpath %.cpp ../../src

SOURCE=bubble.cpp fjapp_ksp.cpp fjapp_bubbleNS.cpp fjlib_cubicspline.cpp fjlib_curve.cpp fjlib_interpolation.cpp fjlib_memini.cpp fjlib_polynomial.cpp fjlib_polyroots.cpp fjlib_smoothing.cpp fjlib_solver.cpp 

#EXTRA_DEP=fjapp_bubble_driver.cpp fjapp_vof.cpp fjapp_vof_staggered.cpp

PROJECT=bubble

OBJECTS=${SOURCE:.cpp=.o}

${PROJECT}: ${OBJECTS} chkopts

	-${CLINKER} -o $@ ${OBJECTS} ${PETSC_KSP_LIB} ${GPROF}

all: ${PROJECT}

