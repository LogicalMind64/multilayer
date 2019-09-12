all: takagi-taupin

C=gcc
CPP=g++
INCS=-I/users/osterhof/.root/include
LIBS=-L/users/osterhof/.root/lib

takagi-taupin: tt-intERGO.c tt-intERGO.h paramsERGO.h
	$(C) -o $@ $< $(INCS) $(LIBS) -O2 -Wall -Wextra -Werror -lm -lneon -lgsl -lgslcblas

.PHONY: clean
clean:
	-rm takagi-taupin plotml

