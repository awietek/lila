CC        = mpicxx
CCARCH    = -cxx=icpc -std=c++11
CCOPT     = -O3 -Wall -Wextra -pedantic
# LAPACKLIB = -llapack -lblas
LAPACKLIB = -mkl -DLILA_USE_MKL
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
