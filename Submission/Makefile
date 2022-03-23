solver: main.cpp
	g++ -o solver main.cpp -lboost_program_options -O3 -O2 -fopenmp

test1: solver
	./solver --Nx 101 --Ny 101 --a 0.75 --b 0.06 --eps 50.0 --mu1 5.0 --mu2 0.0

test2: solver
	./solver --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13.0 --mu1 5.0 --mu2 0.0

test3: solver
	./solver --Nx 101 --Ny 101 --a 0.5 --b 0.1 --eps 50.0 --mu1 5.0 --mu2 0.0 

test4: solver
	./solver --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01

main.o:	main.cpp	
	g++ -o main.o -c main.cpp


.PHONY: clean

clean:
	rm -f *.0 solver
