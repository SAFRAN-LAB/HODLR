CC	=g++
CFLAGS	=-c -Wall -O4 -ffast-math -ffinite-math-only -I header_complex/
LDFLAGS	=
SOURCES	=./examples_complex/HODLR_Test.cpp ./examples_complex/KDTree.cpp
KERNEL	=-DSINC # use -DEXPONENTIAL, -DGAUSSIAN, -DSINC, -DQUADRIC, -DINVERSEQUADRIC, -DMULTIQUADRIC, -DINVERSEMULTIQUADRIC, -DR2LOGR, -DLOGR, -DONEOVERR
DIM	=-DONE  # use -DONE, -DTWO, -DTHREE
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./exec_complex/HODLR_Test_complex

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(DIM) $< -o $@

clean:
	rm -rf *.out ./examples_complex/*.o ./exec_complex/*

tar:
	tar -zcvf HODLR.tar.gz ./makefile_complex.mk ./exec_complex ./header_complex ./examples_complex ./README.md ./LICENSE.md
