CC        = g++
CCARCH    = -std=c++11
CCOPT     = -O3 -Wall -Wextra -pedantic -ferror-limit=1
# LAPACKLIB = -llapack -lblas
# LAPACKLIB = -lmkl_rt -DLILA_USE_MKL
LILAD =  -DLILA_USE_ACCELERATE
LAPACKLIB =  -framework accelerate
INCLUDE   = -I..

EXECUTABLE = tests
OBJECTS    = tests.o test_matrix.o test_special.o test_random.o test_mult.o test_vector.o test_solve.o test_matrixstring.o test_qr.o test_eigen.o test_matrixfunction.o test_tmatrix.o test_lanczos.o test_dot.o test_resize.o test_bandlanczos.o test_gramschmidt.o test_expm.o test_cholesky.o test_lobpcg.o test_expsymv.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CCARCH) $(CCOPT) -o $(EXECUTABLE) $(OBJECTS) $(LAPACKLIB)   

%.o: %.cpp
	$(CC) $(CCARCH) $(CCOPT) $(LILAD) $(INCLUDE) -c $< -o $@ 
clean:
	rm *.o

rebuild: clean $(EXECUTABLE)
