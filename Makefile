
FC = gfortran  
FLAGS = -O2
#FLAGS = -g -CB -check all 
#FLAGS = -g -Warray-bounds -fbounds-check
 
TARGET1 = ns

SOURCES1 = common.f90 \
	   solvers.f90 \
	   main.f90

OBJS1 = $(SOURCES1:.f90=.o)

%.o: %.f90 
	$(FC) $(FLAGS) -c  $<

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(FC) -o $(TARGET1) $(OBJS1) $(LIBS)



clean:
	rm *.o *.mod $(TARGET1) 

cleandat:
	rm *.dat *.txt *.raw 

solvers.o : common.o
main.o : common.o solvers.o
