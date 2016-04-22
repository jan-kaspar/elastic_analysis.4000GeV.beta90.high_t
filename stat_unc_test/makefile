first: fit_test fit_low_t_range

fit_test: fit_test.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall -Wextra fit_test.cc -o fit_test

fit_low_t_range: fit_low_t_range.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall -Wextra fit_low_t_range.cc -o fit_low_t_range
