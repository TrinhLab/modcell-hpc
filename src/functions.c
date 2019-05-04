/* Genetic operators and functions to compute design objectives  */
/* In general methods here should be MOEA independent */

#include "modcell.h"


// operators+fitness.c?
/* Genetic operators*/
/* Objective calculation*/


//WIP: Temporary function
double
calculate_objectives(MCproblem *mcp, Individual *indv){
    for (int k = 0; k < mcp->n_models; k++)
        indv->objectives[k] = 0;
}
