main: module_parmetis.hpp main.cpp makefile
	mpic++ -g -O3 -c module_parmetis.hpp main.cpp -lparmetis -lmetis
	mpic++ -g -O3 module_parmetis.hpp main.cpp -o main -lparmetis -lmetis

clean:
	rm *gnuplot*
