CC = gcc
LIBS = -lglpk
CFLAGS = -O3 -march=native

ifeq ($(flags), develop)
	CFLAGS = -Wall -Wpedantic -g -Wextra -Wno-sign-compare
endif
ifeq ($(flags), optimize)
	CFLAGS = -O3 -march=native -Wall
endif

ifeq ($(link), static)
	LIBS = ../notes/glpk-4.65/src/.libs/libglpk.a -lm
endif
ifeq ($(link), dynamic)
	LIBS = -lglpk
endif

SRC = $(wildcard *.c)
ODIR = obj

_OBJ = $(SRC:.c=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

modcell: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f $(ODIR)/*.o modcell

.PHONY : clean
