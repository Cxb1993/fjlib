CC_LINKER=mpiCC
GPROF=
CFLAGS=
FFLAGS=
CPPFLAGS=-w -DNDEBUG -I../../include
PPPFLAGS=
vpath %.cpp ../../src

include ${PETSC_DIR}/bmake/common/base

SOURCE=sphere.cpp fjapp_SphereAdsorp.cpp fjapp_ksp.cpp fjapp_kspm.cpp fjlib_cubicspline.cpp fjlib_curve.cpp fjlib_interpolation.cpp fjlib_polynomial.cpp fjlib_polyroots.cpp fjlib_solver.cpp fjlib_root.cpp

#EXTRA_DEP=fjapp_bubble_driver.cpp fjapp_vof.cpp fjapp_vof_staggered.cpp

PROJECT=sphere

OBJECTS=${SOURCE:.cpp=.o}

${PROJECT}: ${OBJECTS} chkopts

	-${CLINKER} -o $@ ${OBJECTS} ${PETSC_KSP_LIB}

all: ${PROJECT}

