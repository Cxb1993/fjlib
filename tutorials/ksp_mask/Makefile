CC_LINKER=mpiCC
GPROF=
CFLAGS=
FFLAGS=
CPPFLAGS=-w -DNDEBUG -O2 -I../../include
PPPFLAGS=

include ${PETSC_DIR}/bmake/common/base

SOURCE=kspm.cpp fjapp_ksp.cpp fjapp_kspm.cpp

PROJECT=kspm

OBJECTS=${SOURCE:.cpp=.o}

${PROJECT}: ${OBJECTS} chkopts

	-${CLINKER} -o $@ ${OBJECTS} ${PETSC_KSP_LIB}

all: ${PROJECT}

