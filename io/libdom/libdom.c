#include <stdlib.h>
#include <math.h>
#include "libdom.h"

#define OBJTOL 0.0000001 /* Differences in objectives below this value are ignored. This is because floating point precission cannot be trusted, and not some practical threshold. */

void
fdom(const int m, const int n, const double **x, double *y) // TODO: y should be int
{
    size_t i, ii;

    for(i=0; i < m; i++) /* All are non-dominated until proven otherwise */
        y[i] = 0;

    for(i=0; i < m; i++) { /* find if solution i is dominated */
        for(ii=0; ii < m; ii++) { //Currently comparisons are not comutable (i.e., i vs ii, is not the same as ii vs i, if this property is introduced then we can begin at ii=i+1
            if( (i != ii) && (y[ii] != 1) ) { /* Do not compare against already dominated solutions, since the solution that dominate them can also dominate this one */

                /* If objective values are equal in all cases, the solutions are not comparable */
                if (!is_equal(n, x, i, ii)) {
                    if(is_i_dominated_by_ii(n, x, i, ii)) {
                        y[i] = 1;
                        break; /* Break out of ii loop, since solution i has been proven to be dominated */
                    }
                }

            }
        }
    }
}

int
is_equal(const int n, const double **x, const int i, const int ii)
{
    size_t j;

    for(j=0; j < n; j++) {
        if( fabs(x[i][j] - x[ii][j]) > OBJTOL ) {
            return 0;
        }
    }
    return 1;
}

int
is_i_dominated_by_ii(const int n, const double **x, const int i, const int ii)
{
    /* Assumes solutions are not equal */
    size_t j;

    for(j=0; j < n; j++) {
        if( x[i][j] > x[ii][j] ) { // TODO: Should OBJTOL be applyed here?
            return 0;
         }
     }
    return 1;

}
