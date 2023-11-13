CC=mpif90

CFLAGS=-O3 -std=legacy -W

init: $(wildcard *.for)
	$(LINK.c) -I. \
	$^ $(LOADLIBS) $(LDLIBS) -o $@
	@rm -f *.o

all:  init

clean:
	rm -rf init test.cp
