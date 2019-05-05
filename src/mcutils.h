/* General purpose routines*/

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LMAX 255

const char *
get_file_extension(const char *fspec) {
    char *e = strrchr (fspec, '.');
    if (e == NULL)
        e = "";
    return ++e; /* pointer increased to avoid "." in extension */
}

char *
remove_file_extension(char* mystr) {
    char *retstr;
    char *lastdot;
    if (mystr == NULL)
         return NULL;
    if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
        return NULL;
    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, '.');
    if (lastdot != NULL)
        *lastdot = '\0';
    return retstr;
}

int
is_extension(char *file_name, char *extension){
    if(strcmp(get_file_extension(file_name), extension) == 0)
        return 1;
    return 0;
}

void
set_full_path(char *full_path, const char *parent_path, const char *file_path) {
    	    strcpy(full_path, parent_path);
    	    strcat(full_path, file_path);
}

/* A dynamic array of char* to get file or directory contents*/
typedef struct  {
	char **array;
	size_t n;
} Charlist;

void
alloc_charlist(Charlist *charlist){
if (!(charlist->array = calloc (LMAX, sizeof *charlist->array ))) {
	fprintf (stderr, "error: memory allocation failed.");
	exit(-1);
}
}

size_t
realloc_charlist(Charlist *charlist, size_t *p_lmax){
	size_t lmax = (size_t)&p_lmax;
        if (charlist->n == lmax) {      /* if lmax lines reached, realloc   */
            char **tmp = realloc (charlist->array, lmax * 2 * sizeof *charlist->array);
            if (!tmp)
		exit(-1);
            charlist->array = tmp;
            return lmax * 2;
        }
	return lmax;
}

void free_charlist(Charlist c){
    	for (int i = 0; i < c.n; i++)
        	free (c.array[i]);
    	free (c.array);
}

/* Read each file name of a directory*/
Charlist
read_dir(const char *dir_path){
	Charlist charlist = {NULL, 0};
    	size_t it = 0;
    	size_t lmax = LMAX;         /* current array pointer allocation */
    	DIR *d;
    	struct dirent *dir;

    	d = opendir(dir_path);
    	if (!d) {
        	fprintf (stderr, "error: opening problem directory: '%s'.", dir_path);
		exit(-1);
    	}

	alloc_charlist(&charlist);
    	while ((dir = readdir(d)) != NULL) {
        	charlist.array[charlist.n++] = strdup(dir->d_name);
		lmax = realloc_charlist(&charlist, &lmax);
    	}
    	closedir(d);
    	return charlist;
}

/* NOTE: This method may have issues dealing with large files */
Charlist
read_file(const char *file_path){
	Charlist charlist = {NULL, 0};
	char *ln = NULL;
    	size_t n = 0;               /* buf size, 0 use getline default  */
    	ssize_t nchr = 0;           /* number of chars actually read    */
    	size_t it = 0;
    	size_t lmax = LMAX;         /* current array pointer allocation */
    	FILE *fp = NULL;


	if (!(fp = fopen (file_path, "r"))) {
        	fprintf (stderr, "error: file open failed '%s'.", file_path);
        	exit(-1);
    	}

	alloc_charlist(&charlist);
    	while ((nchr = getline (&ln, &n, fp)) != -1)  /* Note:- ssize_t getline (char **ln, size_t *n, FILE *fp) */
    	{
        	while (nchr > 0 && (ln[nchr-1] == '\n' || ln[nchr-1] == '\r'))
            		ln[--nchr] = 0;     /* strip newline or carriage rtn    */

        	charlist.array[charlist.n++] = strdup (ln);
		lmax = realloc_charlist(&charlist, &lmax);
    	}

    	if (fp) fclose (fp);
    	if (ln) free (ln);

    	return charlist;
}

