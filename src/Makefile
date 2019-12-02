CC=mpicc
LIBS = ../bin/libglpk.a -lm
MODCELL_V_STRING := $(shell git rev-parse HEAD | sed 's:\(.*\):\x27"\1"\x27:')
CFLAGS = -O3 -DMODCELL_V_STRING=$(MODCELL_V_STRING) -Wall


ifeq ($(flags), develop)
	CFLAGS = -Wall -Wpedantic -g -Wextra -Wno-sign-compare -DMODCELL_V_STRING=$(MODCELL_V_STRING)
endif
ifeq ($(flags), optimize)
	CFLAGS = -O3 -march=native -Wall -DMODCELL_V_STRING=$(MODCELL_V_STRING)
endif


ifeq ($(flags), portable)
	CFLAGS = -O3 -Wall -DMODCELL_V_STRING=$(MODCELL_V_STRING)
	LIBS = ../bin/libglpk.a -lm
endif

ifeq ($(link), static)
	LIBS = ../bin/libglpk.a -lm
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
