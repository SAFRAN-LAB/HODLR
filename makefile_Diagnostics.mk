CC	=g++
CFLAGS	=-c -Wall -O4 -ffast-math -ffinite-math-only -I header/
LDFLAGS	=
SOURCES	=./examples/HODLR_diagnostics.cpp ./examples/KDTree.cpp
KERNEL	=-DGAUSSIAN # use -DEXPONENTIAL, -DGAUSSIAN, -DSINC, -DQUADRIC, -DINVERSEQUADRIC, -DMULTIQUADRIC, -DINVERSEMULTIQUADRIC, -DR2LOGR, -DLOGR, -DONEOVERR
DIM	=-DTWO  # use -DONE, -DTWO, -DTHREE
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./exec/HODLR_Diagnostics

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(DIM) $< -o $@

clean:
	rm -rf *.out ./examples/*.o ./src/*.o ./exec/*

tar:
	tar -zcvf HODLR.tar.gz ./makefile_HODLR_Test.mk ./exec ./src ./header ./examples ./README.md ./LICENSE.md
