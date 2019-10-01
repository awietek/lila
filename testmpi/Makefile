CC        = mpicxx
CCARCH    =  -std=c++11
CCOPT     = -O3 -Wall -Wextra -pedantic
# LAPACKLIB = -llapack -lblas
# LAPACKLIB = -lmkl_rt -DLILA_USE_MKL
LAPACKLIB = -framework Accelerate
INCLUDE   = -I..

EXECUTABLE = tests
OBJECTS    = tests.o  test_mpi.o test_dot.o test_add.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CCARCH) $(CCOPT) -o $(EXECUTABLE) $(OBJECTS) $(LAPACKLIB)   

%.o: %.cpp
	$(CC) $(CCARCH) $(CCOPT) $(LAPACKLIB) $(INCLUDE) -c $< -o $@ 
clean:
	rm *.o

rebuild: clean $(EXECUTABLE)
