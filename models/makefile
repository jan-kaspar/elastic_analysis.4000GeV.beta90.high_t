first: sample_models

sample_models: sample_models.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall -Wextra sample_models.cc -o sample_models
