# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

modcell: modcell.o pcg_basic.o
	gcc modcell.o pcg_basic.o -lglpk -lm -o modcell

pcg_basic.o: pcg_basic.c
	gcc -c pcg_basic.c

modcell.o: modcell.c
	gcc -c modcell.c

clean:
	rm -f *.o modcell
