OBJ    = main.o solverwa.o solver.o mq.o io.o mff.o hydro.o eos.o matrix.o
CC=g++
CFLAGS=-W -Wall -Wextra -pedantic -pg -O3 -c -std=c++11
all: exe

exe: $(OBJ)
	$(CC) -fopenmp -pg -o exe $(OBJ)

main.o: main.cpp solver.hpp solverwa.hpp mq.hpp eos.hpp
	$(CC) -fopenmp $(CFLAGS) main.cpp

solver.o: solver.cpp solver.hpp mff.hpp io.hpp
	$(CC) $(CFLAGS) solver.cpp

solverwa.o: solverwa.cpp solverwa.hpp solver.hpp
	$(CC) $(CFLAGS) solverwa.cpp

mq.o: mq.cpp mq.hpp matrix.hpp io.hpp
	$(CC) $(CFLAGS) mq.cpp

io.o: io.cpp io.hpp matrix.hpp
	$(CC) $(CFLAGS) io.cpp

mff.o: mff.cpp mff.hpp hydro.hpp eos.hpp
	$(CC) $(CFLAGS) mff.cpp

hydro.o: hydro.cpp matrix.hpp hydro.hpp eos.hpp
	$(CC) $(CFLAGS) hydro.cpp
eos.o: eos.cpp matrix.hpp
	$(CC) $(CFLAGS) eos.cpp

matrix.o: matrix.cpp matrix.hpp
	$(CC) $(CFLAGS) matrix.cpp

