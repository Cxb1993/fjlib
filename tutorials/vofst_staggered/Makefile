CC=g++
CFLAGS=-w -I../../include

vpath %.cpp ../../src

SOURCE=vofst.cpp fjlib_cubicspline.cpp fjlib_curve.cpp fjlib_interpolation.cpp fjlib_polynomial.cpp fjlib_polyroots.cpp fjlib_solver.cpp 
OBJECTS=${SOURCE:.cpp=.o}

PROJECT=vofst

.PREFIXES= .cpp .o

.cpp.o:
	${CC} -c ${CFLAGS} $< -o $@ -g 

${PROJECT}: ${OBJECTS}
	${CC} -o $@ ${OBJECTS} -lm ${LINKFLAGS}

all: ${PROJECT}

