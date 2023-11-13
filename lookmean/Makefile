# CC=mpif90
CC=gfortran
# CFLAGS=-O3

lookmean: $(wildcard *.for )
	$(LINK.c) -I. \
	$^ $(LOADLIBS) $(LDLIBS) -o $@

all:  lookmean

clean:
	rm -rf lookmean
