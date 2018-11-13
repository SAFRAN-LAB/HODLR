CC		   =g++-6
CFLAGS	   =-c -I $(EIGEN_PATH) -fopenmp -Wall -Ofast -funroll-loops -ffast-math -ffinite-math-only -I src/
LDFLAGS	   =-fopenmp -I $(EIGEN_PATH)
SOURCES	   =./examples/testHODLR.cpp
OBJECTS	   =$(SOURCES:.cpp=.o)
EXECUTABLE =./exec/HODLR_Test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.out ./examples/*.o ./exec/*
