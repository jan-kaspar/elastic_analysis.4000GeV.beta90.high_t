all: simu make_generators

simu: simu.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall -Wextra simu.cc -o simu

make_generators: make_generators.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall -Wextra make_generators.cc -o make_generators
