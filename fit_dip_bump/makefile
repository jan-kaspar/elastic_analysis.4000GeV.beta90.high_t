all: make_fits

make_fits: make_fits.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall -Wextra -g\
		make_fits.cc -o make_fits
