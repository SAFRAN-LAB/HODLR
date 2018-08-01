CC		=g++
CFLAGS		=-c -Wall -Ofast -fopenmp -ffast-math -ffinite-math-only
LDFLAGS		=-Ofast -fopenmp
SOURCES		=./testHODLR.cpp ./HODLR.hpp
KERNEL		=-DLOGR	# use -DINVEXPR, -DLOGR, -DONEOVERR
OBJECTS		=$(SOURCES:.cpp=.o)
EXECUTABLE	=./testHODLR

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $(KERNEL) $(HOMOG) $< -o $@

clean:
	rm a.out testHODLR *.o
