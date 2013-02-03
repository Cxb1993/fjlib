CC=g++
CFLAGS=-w -I../../include 
#LINKFLAGS=-lboost_unit_test_framework

SRC1=parser.cpp
SRC2=chart.cpp

vpath %.cpp ../../src
vpath %.h ../../include

PROJECT=parser

OBJ1=${SRC1:.cpp=.o}
OBJ2=${SRC2:.cpp=.o}

.PREFIXES= .cpp .o

.cpp.o: 
	${CC} -c ${CFLAGS} $< -o $@ -g

parser: ${OBJ1}
	${CC} -o $@ ${OBJ1} -lm ${LINKFLAGS} -g

chart: ${OBJ2}
	${CC} -o $@ ${OBJ2} -lm ${LINKFLAGS} -g

