# Makefile created with Create Makefile.sh 16/02/2012

COPT      = -O3
CFLAGS    = $(COPT) -ftree-vectorize -funroll-loops -fopenmp -Wall -Wextra -ansi
CSTRICT   = -Wall -Wextra -ansi
CXXOPT    = -O3
CXXFLAGS  = $(CXXOPT) -ftree-vectorize -funroll-loops -fopenmp -Wall -Wextra -ansi
CXXSTRICT = -Wall -Wextra -ansi
LDFLAGS   = -lpng -lgomp
EXEC      = ksvd

default: $(EXEC)
all: $(EXEC)

# ------- C files -------
./mt19937ar.o: ./mt19937ar.c ./mt19937ar.h
	$(CC) $(CFLAGS) -c ./mt19937ar.c -o ./mt19937ar.o

./io_png.o: ./io_png.c ./io_png.h
	$(CC) $(CFLAGS) -c ./io_png.c -o ./io_png.o

# ------- C++ files -------
./main.o: ./main.cpp
	$(CXX) $(CXXFLAGS) -c ./main.cpp -o ./main.o

./lib_svd.o: ./lib_svd.cpp ./lib_svd.h
	$(CXX) $(CXXFLAGS) -c ./lib_svd.cpp -o ./lib_svd.o

./addnoise_function.o: ./addnoise_function.cpp ./addnoise_function.h
	$(CXX) $(CXXFLAGS) -c ./addnoise_function.cpp -o ./addnoise_function.o

./lib_ormp.o: ./lib_ormp.cpp ./lib_ormp.h
	$(CXX) $(CXXFLAGS) -c ./lib_ormp.cpp -o ./lib_ormp.o

./utilities.o: ./utilities.cpp ./utilities.h
	$(CXX) $(CXXFLAGS) -c ./utilities.cpp -o ./utilities.o

./ksvd.o: ./ksvd.cpp ./ksvd.h
	$(CXX) $(CXXFLAGS) -c ./ksvd.cpp -o ./ksvd.o

# ------- Main -------
$(EXEC):  ./mt19937ar.o ./io_png.o ./main.o ./lib_svd.o ./addnoise_function.o ./lib_ormp.o ./utilities.o ./ksvd.o
	$(CXX)  ./mt19937ar.o ./io_png.o ./main.o ./lib_svd.o ./addnoise_function.o ./lib_ormp.o ./utilities.o ./ksvd.o $(LDFLAGS) -o $(EXEC)

lint: 
	$(MAKE) CFLAGS="$(CFLAGS) $(CSTRICT)" CXXFLAGS="$(CXXFLAGS) $(CXXSTRICT)"

clean: 
	rm -f *.o

distclean: clean
	rm -f $(EXEC)

