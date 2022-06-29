objects = main.o peblock.o  pecsb2.o Triangle.o
pe:$(objects)
	g++ -O4 -mavx2 -fopenmp -o pe $(objects) -I /usr/local/hdf5/ -L /usr/local/hdf5/lib /usr/local/hdf5/lib/libhdf5.a -lsz -lz -lm -ldl
main.o:pecsb2.h
	g++ -O4 -mavx2 -fopenmp -o main.o -c main.cpp
peblock.o:peblock.h DyArray.h Solve_Tridiagonal.h
	g++ -O4 -mavx2 -fopenmp -o peblock.o -c peblock.cpp
pehor2.o:pehor2.h DyArray.h Solve_Tridiagonal.h
	g++ -O4 -mavx2 -fopenmp -o  
pecsb2.o:DyArray.h  peblock.h Triangle.h
	g++ -fopenmp -o pecsb2.o -c pecsb2.cpp
Triangle.o:Triangle.h
	g++ -fopenmp -o Triangle.o -c Triangle.cpp
clean:
	rm $(objects)

