#GPROF=-g -pg
CC_LINKER=mpiCC
GPROF=
CFLAGS=
FFLAGS=
CPPFLAGS=-w -DNDEBUG -O2 -I../../include ${GPROF}
PPPFLAGS=

include ${PETSC_DIR}/bmake/common/base

SOURCE=bubble.cpp fjapp_ksp.cpp fjapp_bubbleNS.cpp fjapp_MInterface.cpp fjapp_bubbleNS_MI.cpp fjlib_cubicspline.cpp fjlib_curve.cpp fjlib_interpolation.cpp fjlib_memini.cpp fjlib_polynomial.cpp fjlib_polyroots.cpp fjlib_smoothing.cpp fjlib_solver.cpp 

#EXTRA_DEP=fjapp_bubble_driver.cpp fjapp_vof.cpp fjapp_vof_staggered.cpp

PROJECT=bubble

OBJECTS=${SOURCE:.cpp=.o}

${PROJECT}: ${OBJECTS} chkopts

	-${CLINKER} -o $@ ${OBJECTS} ${PETSC_KSP_LIB} ${GPROF}

all: ${PROJECT}

