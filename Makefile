CC=gcc
OUTPUT_OPTION=-Wall -O2 -ferror-limit=50 -MMD -MP -o $@
LDLIBS=-lm
 
SRC=$(wildcard src/*.c)
OBJ=$(SRC:.c=.o)
DEP=$(SRC:.c=.d)
 
.PHONY: clean
 
sumatra: $(OBJ)
	$(LINK.o) -O2 $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	rm -f $(OBJ) $(DEP) sumatra

-include $(DEP)
