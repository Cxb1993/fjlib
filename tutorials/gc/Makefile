CC=g++
CFLAGS=-w -I../../include 
LINKFLAGS=-lboost_unit_test_framework

SOURCE=gc.cpp fjlib_cubicspline.cpp fjlib_curve.cpp fjlib_interpolation.cpp fjlib_polynomial.cpp fjlib_polyroots.cpp fjlib_solver.cpp
EXTDEP=fjapp_gc.h fjapp_vof.h

vpath %.cpp ../../src
vpath %.h ../../include

PROJECT=gc

OBJECTS=${SOURCE:.cpp=.o}


.PREFIXES= .cpp .o

.cpp.o: ${EXTDEP}
	${CC} -c ${CFLAGS} $< -o $@ -g

${PROJECT}: ${OBJECTS} ${EXTDEP}
	${CC} -o $@ ${OBJECTS} -lm ${LINKFLAGS} -g

all: ${PROJECT}

