CC=mpif90

CFLAGS=-O3 -std=legacy -W

PipeModified: $(wildcard *.for)
	$(LINK.c) -I. \
	$^ $(LOADLIBS) $(LDLIBS) -o $@
	@rm -f *.o

all:  PipeModified

clean:
	rm -rf PipeModified
