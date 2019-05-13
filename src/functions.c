/* The methods in this file MOEA independent abstractions (although some variables that appear in the datastructures specific to NSGA-II may also appear here). That includes:
 * Genetic operators, functions to compute design objectives, and other MOEA functions such as checking domination, etc.
 */

#include <stdlib.h>
#include "modcell.h"

extern glp_smcp param;

void copy_individual(MCproblem *mcp, Individual *indv_source, Individual *indv_dest);
int find_domination(MCproblem *mcp, Individual *indv_a, Individual *indv_b);
void crossover(MCproblem *mcp, Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
void mutation(MCproblem *mcp, Individual *indv);
void enforce_module_constraints(MCproblem *mcp, Individual *indv);
void calculate_objectives(MCproblem *mcp, Individual *indv);


void
copy_individual(MCproblem *mcp, Individual *indv_source, Individual *indv_dest)
{
    int i, j, k;

    for (j=0; j < mcp->n_vars; j++)
        indv_dest->deletions[j] = indv_source->deletions[j];

    for (k=0; j < mcp->n_models; k++)
        for (j=0; j < mcp->n_vars; j++)
            indv_dest->modules[k][j] = indv_source->modules[k][j];

    for (i=0; i < mcp->n_models; i++) {
        indv_dest->objectives[i] = indv_source->objectives[i];
        indv_dest->penalty_objectives[i] = indv_source->penalty_objectives[i];
    }
    indv_dest->rank = indv_source->rank;
    indv_dest->crowding_distance = indv_source->crowding_distance;
}

/* Checks if indv1 dominates indv2. Assumes objectives are maximized! Returns:
   1DOMINATES2 if indv1 dominates indv2
   2DOMINATES1 if indv2 dominates indv1
   NONDOMINATED if both are non dominated
*/
int
find_domination(MCproblem *mcp, Individual *indv_a, Individual *indv_b)
{
    int a_dominates_b = 1;
    int b_dominates_a = 1;

    for (int i=0; i < mcp->n_models; i++){
        if(indv_a->penalty_objectives[i] > indv_b->penalty_objectives[i])
            b_dominates_a = 0; // can stop if a_dominates_b==0
        if(indv_b->penalty_objectives[i] > indv_a->penalty_objectives[i])
            a_dominates_b = 0; //can stop if 2dominates1==0
    }
    if(a_dominates_b & b_dominates_a) return NONDOMINATED; /* both objective vectors are equal */
    if(a_dominates_b) return A_DOMINATES_B;
    if(b_dominates_a) return B_DOMINATES_A;
    return NONDOMINATED;
}

/* Two point binary crossover of two individuals
   - The crossover probability  is evaluated here and if crossoverr is not perform the childs will match the parents
   - Crossover on module reactions is done on each model indepently. However, the  crossover  sites are the same that in deletions, given the relation between both variables this is a better way to preserve blocks. This is tricky since it might also be good to be able to get rid of modules.
   -TODO: Could refactor for contiguous memory access (e.g.  for [k][j] variables, loop through j within loop through k) but this would lead to additional for loops. Although currently if(mcp->beta>0) is checked redundantly.
*/

#define  FILL \
    child1->deletions[j] = parent1->deletions[j]; \
    child2->deletions[j] = parent2->deletions[j];

#define  FILLB \
    if (mcp->beta > 0) { \
        for (k=0; k < mcp->n_models; k++) { \
            child1->modules[k][j] = parent1->modules[k][j]; \
            child2->modules[k][j] = parent2->modules[k][j]; \
        }\
    }

void
crossover(MCproblem *mcp, Individual *parent1, Individual *parent2, Individual *child1, Individual *child2)
{
    int j, k, temp, site1, site2;

    if ( (double)pcg32_boundedrand(100)/100 <= mcp->crossover_probability)  {
        site1 = pcg32_boundedrand(mcp->n_vars);
        site2 = pcg32_boundedrand(mcp->n_vars);
        if (site1 > site2) { /* swap variables */
            temp = site1;
            site1 = site2;
            site2 = temp;
        }
        for (j=0; j < site1; j++) {
            FILL
            FILLB
        }
        for (j=site1; j < site2; j++){
            child1->deletions[j] = parent2->deletions[j];
            child2->deletions[j] = parent1->deletions[j];
            if (mcp->beta > 0) {
                for (k=0; k < mcp->n_models; k++) {
                    child1->modules[k][j] = parent2->modules[k][j];
                    child2->modules[k][j] = parent1->modules[k][j];
                }
            }
        }
        for (j=site2; j < mcp->n_vars; j++){
            FILL
            FILLB
        }
    } else { /* No crossover is done */
        for (j=0; j < mcp->n_vars; j++) {
            FILL
            FILLB
        }
    }
}


/* Binary mutation of individual
- A random bit might be flipped in deletion array and each module reaction array independently. Flipping the same bit for deletions and all modules would be useless.
*/
void
mutation(MCproblem *mcp, Individual *indv)
{
    int k, site;

    if ( (double)pcg32_boundedrand(100)/100 <= mcp->mutation_probability)  {
        site = pcg32_boundedrand(mcp->n_vars);
        indv->deletions[site] = !indv->deletions[site];
    }

    if (mcp->beta > 0) {
        for (k=0; k < mcp->n_models; k++) {
            if ( (double)pcg32_boundedrand(100)/100 <= mcp->mutation_probability)  {
                site = pcg32_boundedrand(mcp->n_vars);
                indv->modules[k][site] = !indv->modules[k][site];
            }
        }
    }
}


/* After crossover are mutation are done, they may generate individuals that violate the two module reaction related constraints. This method enforces both constraints as follows:
   1. Removes modules that are not deletions
   2. If number of modules is above limit (beta), randomly removes modules until within limit.
Notes:
    - Refactor with linked lists?
    - pcg32 does not provide a method to obtain a list of non-repeated random numbers.
*/
void
enforce_module_constraints(MCproblem *mcp, Individual *indv)
{
    int i, j, k, n_module_rxn, module_diff, n_removed_module, target;
    int module_rxn_idx[MAX_MODULES] = {-1}, is_removed_module[MAX_MODULES] = {-1};

   /*  Removes modules that are not deletions */
    for (k=0; k < mcp->n_models; k++) {
        for (j=0; j < mcp->n_vars; j++) {
            if ( (indv->deletions[j] == !DELETED_RXN) & (indv->modules[k][j] == MODULE_RXN))
                indv->modules[k][j] = !MODULE_RXN;
        }
    }

    for (k=0; k < mcp->n_models; k++) {
        /* Identify module reactions */
        n_module_rxn = 0;
        for (j=0; j < mcp->n_vars; j++) {
            if (indv->modules[k][j] == MODULE_RXN) {
                module_rxn_idx[n_module_rxn] = j;
                n_module_rxn++;
            }
        }
        /* Randomly remove additional modules  */
        module_diff = n_module_rxn - mcp->beta;
        if(module_diff > 0){
            for (i=0; i < module_diff; i++) /* Reset removed modules */
                is_removed_module[i] = 0;
            n_removed_module = 0;
            while (module_diff != n_removed_module){
                target = pcg32_boundedrand(module_diff);
                if (is_removed_module[target] == 0){
                    is_removed_module[target] = 1;
                    n_removed_module++;
                    indv->modules[k][module_rxn_idx[target]] = 0; /* apply removal */
                }
            }
        }
    }
}

/*
   This function will set indv->objectives and indv->penalty_objectives

   Notes:
   - Currently the objective type is specified as part of the input file, this method does not do any manipulation of the models to set the appropriate design objective.
   - Currently only computes wGCP.

    Optimizations (perform them prior to calling this function?):
        - Look up hash table of known solutions
*/
void
calculate_objectives(MCproblem *mcp, Individual *indv)
{
    LPproblem *lp;
    int j,k,n_deletions=0;
    int *change_bound = malloc(mcp->n_vars * sizeof(int));

    /* Preliminary evaluation */
    for (j=0; j < mcp->n_vars; j++)
        if(indv->deletions[j] == 0)
            n_deletions++;

    if (n_deletions == 0){ /* Avoid further evaluation */
        for (k=0; k < mcp->n_models; k++) {
            lp = &(mcp->lps[k]);
            indv->objectives[k] = lp->no_deletion_objective;
            indv->penalty_objectives[k] = lp->no_deletion_objective;
        }
        return;
    }

    /* Objective calculation */
    for (k=0; k < mcp->n_models; k++) {
        lp = &(mcp->lps[k]);

        /* Determine what bounds to change */
        for (j=0; j < mcp->n_vars; j++) {
            change_bound[j] = 0;
            if(lp->cand_col_idx[j] == NOT_CANDIDATE)
                continue;
            if(indv->deletions[j] == 0) {
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

