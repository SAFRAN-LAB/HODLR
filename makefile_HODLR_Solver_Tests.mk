CC	=g++
CFLAGS	=-c -Wall -O3 -ffinite-math-only -I ~/Dropbox/Eigen/
LDFLAGS	=
SOURCES	=HODLR_Solver.cpp HODLR_Solver_Tests.cpp
KERNEL	=-DSINC  # use -DONEOVERR, -DEXPONENTIAL, -DLOGARITHM, -DGAUSSIAN, or -DSINC
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE=HODLR_Solver_Tests

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(KERNEL) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $< -o $@

clean:
	rm -rf *.out HODLR_Solver_Tests *.o

tar:
	tar -zcvf *.cpp *.hpp makefile_HODLR_Solver.mk makefile_HODLR_Solver_Tests.mk
