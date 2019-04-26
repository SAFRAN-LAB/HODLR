CC	        =g++
FC 	        =gfortran
CFLAGS	    =-c -O2 -ffast-math -ffinite-math-only -I ../../include/ -I $(EIGEN_PATH)
FFLAGS      =-c -O3
LDFLAGS	    =-lstdc++
CSOURCES	=../../src/*.cpp
FSOURCES 	=example.f
COBJECTS	=$(CSOURCES:.cpp=.o)
FOBJECTS	=$(FSOURCES:.f=.o)
EXECUTABLE	=example
all: $(FSOURCES) $(CSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(FOBJECTS) $(COBJECTS)
	$(FC) $(FOBJECTS) $(COBJECTS) -o $@ $(LDFLAGS) 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.f.o:
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -rf $(COBJECTS) $(FOBJECTS) $(EXECUTABLE)
