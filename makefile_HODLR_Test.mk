CC	=g++
CFLAGS	=-c -Wall -O3 -ffinite-math-only -I ~/Dropbox/Eigen/
LDFLAGS	=
SOURCES	=KDTree.cpp get_Matrix.cpp partial_Piv_LU.cpp HODLR_Test.cpp HODLR_Node.cpp HODLR_Tree.cpp
KERNEL	=-DSINC  # use -DEXPONENTIAL, -DGAUSSIAN, -DSINC, -DQUADRIC, -DINVERSEQUADRIC, -DMULTIQUADRIC, -DINVERSEMULTIQUADRIC
DIM	=-DTWO  # use -DONE, -DTWO, -DTHREE
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=HODLR_Test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(DIM) $< -o $@

clean:
	rm -rf *.out HODLR_Test *.o *.*~

tar:
	tar -zcvf *.cpp *.hpp makefile_HODLR_Test.mk
