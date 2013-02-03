CC=g++
CFLAGS=-w -I../../include 
LINKFLAGS=-lboost_unit_test_framework

SOURCE=vof_uvp.cpp fjlib_cubicspline.cpp fjlib_curve.cpp fjlib_interpolation.cpp fjlib_polynomial.cpp fjlib_polyroots.cpp fjlib_solver.cpp

PROJECT=vof_uvp

vpath %.cpp ../../src

OBJECTS=${SOURCE:.cpp=.o}

.PREFIXES= .cpp .o

.cpp.o:
	${CC} -c ${CFLAGS} $< -o $@ -g 

${PROJECT}: ${OBJECTS}
	${CC} -o $@ ${OBJECTS} -lm ${LINKFLAGS} -g

all: ${PROJECT}

