GSL_DIR = /afs/cern.ch/cms/slc6_amd64_gcc481/external/gsl/1.10/

INCS = -I$(GSL_DIR)/include
LIBS = -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm

first: \
	numerical_integration\
	unfolding_test

numerical_integration: numerical_integration.cc numerical_integration_common.h common.h
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall numerical_integration.cc -o numerical_integration

unfolding_test: unfolding_test.cc ../../NumericalIntegration.h
	g++ `root-config --libs` -lMinuit `root-config --cflags` $(INCS) $(LIBS) -O3 -Wall unfolding_test.cc -o unfolding_test
