all: make_fits tabulate

make_fits : make_fits.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` -O3 -g \
		make_fits.cc -o make_fits

tabulate : tabulate.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` -O3 -g \
		tabulate.cc -o tabulate
