#

modcell: modcell.o
	gcc modcell.o -lglpk -lm -o modcell

modcell.o: modcell.c
	gcc -c modcell.c
clean:
	rm -f *.o solve
