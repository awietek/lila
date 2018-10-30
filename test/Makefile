CC        = g++
CCARCH    = -std=c++11
CCOPT     = -O3 -Wall -Wextra -pedantic
LAPACKLIB = -llapack -lblas
#LAPACKLIB = -mkl -DPLATFORM_MKL
INCLUDE   = -I..

EXECUTABLE = tests
OBJECTS    = tests.o test_matrix.o test_special.o test_random.o test_mult.o test_vector.o test_solve.o test_matrixstring.o test_qr.o test_eigen.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CCARCH) $(CCOPT) -o $(EXECUTABLE) $(OBJECTS) $(LAPACKLIB)   

%.o: %.cpp
	$(CC) $(CCARCH) $(CCOPT) $(LAPACKLIB) $(INCLUDE) -c $< -o $@ 
clean:
	rm *.o

rebuild: clean $(EXECUTABLE)
