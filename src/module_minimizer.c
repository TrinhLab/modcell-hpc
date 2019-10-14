/* Removes futile module reactions that meet beta constraint but do not enhance the design objective
 * Interface is the same as main moea.c method
 */

#include <stdlib.h>
#include "modcell.h"

void minimize_mr(MCproblem *mcp, Population *parent_population);
void minimize_mr_indv(MCproblem *mcp, Individual *indv, Individual *tindv);

/* Main method
 */
void
minimize_mr(MCproblem *mcp, Population *parent_population)
{
    Population *test_population = malloc(sizeof(Population));
    allocate_population(mcp, test_population, mcp->population_size);

    for (int i=0; i < mcp->population_size; i++) {
        calculate_objectives(mcp, &(parent_population->indv[i])); /* Determine objectives for individual, since these are not parsed from .pop file */
        copy_individual(mcp, &(parent_population->indv[i]), &(test_population->indv[i]));
        minimize_mr_indv(mcp, &(parent_population->indv[i]), &(test_population->indv[i]));
    }

    free_population(mcp, test_population);
}

/* Minimize modules of an individual
 * Notes:
 *      - This method ignores penalty function for alpha.
 *      - TODO: Add verbose (e.g., indv 1: original mr (sum of all modules) :  new  mr:)
 */
void
minimize_mr_indv(MCproblem *mcp, Individual *indv, Individual *tindv)
{

    int j,k;
    int *change_bound = malloc(mcp->n_vars * sizeof(int));


    /* Minimize modules for each objective independently */
    for (k=0; k < mcp->n_models; k++) {
        /* Go through each module and drop it, if objective deteriorates add it back and continue, otherwise keep it as dropped */
        for (j=0; j < mcp->n_vars; j++) {
            if (indv->modules[k*mcp->n_vars + j] == 1) {
                tindv->modules[k*mcp->n_vars + j] = 0;
                calculate_objective(mcp, tindv, k, change_bound);
                if(tindv->objectives[k] + OBJ_TOL < indv->objectives[k])
                    tindv->modules[k*mcp->n_vars + j] = 1;
            }
        }
        /* Update final individual */
        for (j=0; j < mcp->n_vars; j++) {
            indv->modules[k*mcp->n_vars + j]  = tindv->modules[k*mcp->n_vars + j];
        indv->objectives[k] = tindv->objectives[k];
        indv->penalty_objectives[k] = tindv->objectives[k]; /* Note penalty function is ignored */
        }
    }
    free(change_bound);
}
