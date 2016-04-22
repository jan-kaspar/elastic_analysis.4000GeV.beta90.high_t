all: analyze

analyze: analyze.cc analyze_common.h ../../stat.h
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall analyze.cc -o analyze
