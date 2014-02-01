CC	=g++
CFLAGS	=-c -Wall -DNDEBUG -O4 -ffast-math -ffinite-math-only -I header/
LDFLAGS	=
SOURCES	=./src/KDTree.cpp ./src/get_Matrix.cpp ./src/partial_Piv_LU.cpp ./src/HODLR_Node.cpp ./src/HODLR_Tree.cpp ./examples/HODLR_Test.cpp
KERNEL	=-DGAUSSIAN  # use -DEXPONENTIAL, -DGAUSSIAN, -DSINC, -DQUADRIC, -DINVERSEQUADRIC, -DMULTIQUADRIC, -DINVERSEMULTIQUADRIC, -DR2LOGR, -DLOGR, -DONEOVERR
DIM	=-DONE  # use -DONE, -DTWO, -DTHREE
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./exec/HODLR_Test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(DIM) $< -o $@

clean:
	rm -rf *.out ./examples/*.o ./src/*.o ./exec/*

tar:
	tar -zcvf HODLR.tar.gz ./makefile_HODLR_Test.mk ./exec ./src ./header ./examples ./README.md ./LICENSE.md