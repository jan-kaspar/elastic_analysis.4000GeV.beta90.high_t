all: evaluate_stat_significance

evaluate_stat_significance: evaluate_stat_significance.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall -Wextra evaluate_stat_significance.cc -o evaluate_stat_significance
