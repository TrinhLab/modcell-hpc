/* General purpose functions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LMAX 255

/* Get the extension of a file */
const char *getExt (const char *fspec) {
    char *e = strrchr (fspec, '.');
    if (e == NULL)
        e = ""; // fast method, could also use &(fspec[strlen(fspec)]).
    return e;
}

/* Read each line of a text file of arbitrary size into an array*/
struct file {
	char **array; /* array of pointers to char        */
	size_t nlines;
};
struct file readFile(const char *filepath){
	struct file myfile = {NULL, 0};
	char *ln = NULL;            /* NULL forces getline to allocate  */
    size_t n = 0;               /* buf size, 0 use getline default  */
    ssize_t nchr = 0;           /* number of chars actually read    */
    size_t it = 0;              /* general iterator variable        */
    size_t lmax = LMAX;         /* current array pointer allocation */
    FILE *fp = NULL;            /* file pointer                     */


    if (!(fp = fopen (filepath, "r"))) { /* open file for reading    */
        fprintf (stderr, "error: file open failed '%s'.", filepath);
	exit(-1);
    }

    /* allocate LMAX pointers and set to NULL. Each of the 255 pointers will
       point to (hold the address of) the beginning of each string read from
       the file below. This will allow access to each string with array[x].
    */
    if (!(myfile.array = calloc (LMAX, sizeof *myfile.array ))) {
        fprintf (stderr, "error: memory allocation failed.");
	exit(-1);
    }

    /* prototype - ssize_t getline (char **ln, size_t *n, FILE *fp)
       above we declared: char *ln and size_t n. Why don't they match? Simple,
       we will be passing the address of each to getline, so we simply precede
       the variable with the urinary '&' which forces an addition level of
       dereference making char* char** and size_t size_t *. Now the arguments
       match the prototype.
    */
    while ((nchr = getline (&ln, &n, fp)) != -1)    /* read line    */
    {
        while (nchr > 0 && (ln[nchr-1] == '\n' || ln[nchr-1] == '\r'))
            ln[--nchr] = 0;     /* strip newline or carriage rtn    */

        /* allocate & copy ln to array - this will create a block of memory
           to hold each character in ln and copy the characters in ln to that
           memory address. The address will then be stored in array[nlines].
           (nlines++ just increases nlines by 1 so it is ready for the next address)
           There is a lot going on in that simple: array[nlines++] = strdup (ln);
        */
        myfile.array[myfile.nlines++] = strdup (ln);

        if (myfile.nlines == lmax) {      /* if lmax lines reached, realloc   */
            char **tmp = realloc (myfile.array, lmax * 2 * sizeof *myfile.array);
            if (!tmp)
		exit(-1);
            myfile.array = tmp;
            lmax *= 2;
        }
    }

    if (fp) fclose (fp);        /* close file */
    if (ln) free (ln);          /* free memory allocated to ln  */

    return myfile;
}
void freeFile(struct file f){
    for (int it = 0; it < f.nlines; it++)        /* free array memory    */
        free (f.array[it]);
    free (f.array);
}
