/* Genetic operators and functions to compute design objectives */
/* In general methods here should be MOEA independent */

#include <stdlib.h>
#include "modcell.h"

extern glp_smcp param;



/*
   This function will set indv->objectives and indv->penalty_objectives

   Notes:
   - Currently the objective type is specified as part of the input file, this method does not do any manipulation of the models to set the appropriate design objective.
   - Currently only computes wGCP.

    Optimizations (perform them prior to calling this function?):
        - Look up hash table of known solutions
*/
void
calculate_objectives(MCproblem *mcp, Individual *indv){
    LPproblem *lp;
    int j,k,n_deletions=0;
    int *change_bound = malloc(mcp->n_vars * sizeof(int));

    /* Preliminary evaluation */
    for (j=0; j < mcp->n_vars; j++)
        if(indv->deletions[j] == 0)
            n_deletions++;

    if (n_deletions == 0){ /* Avoid further evaluation */
        for (k=0; k < mcp->n_models; k++){
            lp = &(mcp->lps[k]);
            indv->objectives[k] = lp->no_deletion_objective;
            indv->penalty_objectives[k] = lp->no_deletion_objective;
        }
        return;
    }

    /* Objective calculation */
    for (k=0; k < mcp->n_models; k++){
        lp = &(mcp->lps[k]);

        /* Determine what bounds to change */
        for (j=0; j < mcp->n_vars; j++){
            change_bound[j] = 0;
            if(lp->cand_col_idx[j] == NOT_CANDIDATE)
                continue;
            if(indv->deletions[j] == 0){
                change_bound[j] = 1; /* Reaction deleted in the chassis */
                if(mcp->beta > 0)
                    if (indv->modules[k][j] == 1)
                        change_bound[j] = 0; /* Reaction inserted back as module */
            }
        }

        /* Block bounds */
        for (j=0; j < mcp->n_vars; j++)
            if(change_bound[j])
	        glp_set_col_bnds(lp->P, lp->cand_col_idx[j], GLP_FX, 0, 0);

        /* Calculate objectives */
        glp_simplex(lp->P, &param);
        indv->objectives[k] = glp_get_col_prim(lp->P, lp->prod_col_idx)/lp->max_prod_growth;

        /* Reset bounds */
        for (j=0; j < mcp->n_vars; j++)
            if(change_bound[j])
	        glp_set_col_bnds(lp->P, lp->cand_col_idx[j], lp->cand_col_type[j], lp->cand_og_lb[j], lp->cand_og_ub[j]);

        /* Calculate penalty objectives (note that module reaction constraints are strictly enforced by genetic operators)*/
        if (n_deletions > mcp->alpha)
            indv->penalty_objectives[k] = indv->objectives[k]/n_deletions;
        else
            indv->penalty_objectives[k] = indv->objectives[k];
    }
    free(change_bound);
}

