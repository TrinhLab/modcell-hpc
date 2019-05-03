# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

modcell: modcell.o pcg_basic.o allocation.o
	gcc modcell.o pcg_basic.o allocation.o -lglpk -lm -o modcell

pcg_basic.o: pcg_basic.c
	gcc -c pcg_basic.c

modcell.o: modcell.c
	gcc -c modcell.c

allocate.o: allocation.c
	gcc -c allocation.c
clean:
	rm -f *.o modcell
