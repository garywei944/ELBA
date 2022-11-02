all: find_chimera

find_chimera: find_chimera.cpp kppseq.h
	g++ -O2 -o find_chimera find_chimera.cpp -lz

clean:
	rm find_chimera *.dSYM *.o
