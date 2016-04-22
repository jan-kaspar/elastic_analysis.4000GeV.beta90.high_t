all: test_high_thx

test_high_thx: test_high_thx.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall -Wextra -Werror test_high_thx.cc -o test_high_thx
