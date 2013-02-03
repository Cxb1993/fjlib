CC_LINKER=mpiCC
#GPROF=-g -w 
GPROF=-w -DNDEBUG -O2
CFLAGS=
FFLAGS=
CPPFLAGS=-I../../include ${GPROF}
PPPFLAGS=
vpath %.cpp ../../src

include ${PETSC_DIR}/bmake/common/base

SOURCE=bubbleD.cpp fjapp_bubbleD.cpp fjapp_interface.cpp fjapp_ksp.cpp fjapp_kspm.cpp fjlib_cubicspline.cpp fjlib_curve.cpp fjlib_interpolation.cpp fjlib_smoothing.cpp fjlib_polynomial.cpp fjlib_polyroots.cpp fjlib_solver.cpp fjlib_root.cpp fjlib_memini.cpp

#EXTRA_DEP=fjapp_bubble_driver.cpp fjapp_vof.cpp fjapp_vof_staggered.cpp

PROJECT=bubbleD

vpath %.cpp ../../src

OBJECTS=${SOURCE:.cpp=.o}

${PROJECT}: ${OBJECTS} chkopts

	-${CLINKER} -o $@ ${OBJECTS} ${PETSC_KSP_LIB} ${GPROF}

all: ${PROJECT}

