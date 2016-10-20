CC		=g++-5
CFLAGS	=-c -fopenmp -Wall -Ofast -funroll-loops -ffast-math -ffinite-math-only -I header/
LDFLAGS	=-fopenmp
SOURCES	=./examples/testHODLRSymmetric.cpp
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./exec/testHODLRSymmetric

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.out ./examples/*.o ./src/*.o ./exec/*

tar:
	tar -zcvf HODLR.tar.gz ./makefile.mk ./exec ./src ./header ./examples ./README.md ./LICENSE.md
