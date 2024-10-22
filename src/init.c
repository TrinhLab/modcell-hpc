/* Routines to allocate/free  and initialize objects */

#include <stdlib.h>
#include "modcell.h"

void allocate_population(MCproblem *mcp,  Population *pop, size_t pop_size);
void allocate_individual(MCproblem *mcp,  Individual *indv);
void free_population(MCproblem *mcp, Population *pop);
void free_individual(MCproblem *mcp, Individual *indv);

void set_random_population(MCproblem *mcp, Population *pop);
void set_random_individual(MCproblem *mcp,  Individual *indv);
void set_blank_population(MCproblem *mcp, Population *pop);
void set_blank_individual(MCproblem *mcp,  Individual *indv);

void
allocate_MCproblem(MCproblem *mcp, unsigned int n_models, size_t n_vars)
{
    mcp->n_models = n_models;
    mcp->n_vars = n_vars;

    mcp->individual2id = malloc(n_vars * sizeof *mcp->individual2id);
    mcp->model_names = malloc(n_models * sizeof *mcp->model_names);

    LPproblem *lp;
    mcp->lps = malloc(n_models * sizeof(LPproblem));
    for (int k=0; k < n_models; k++) {
        lp = &(mcp->lps[k]);
        lp->P = glp_create_prob();
        lp->cand_col_idx = malloc(n_vars * sizeof(*lp->cand_col_idx));
        lp->cand_og_lb = malloc(n_vars * sizeof(*lp->cand_og_lb));
        lp->cand_og_ub = malloc(n_vars * sizeof(*lp->cand_og_ub));
        lp->cand_col_type = malloc(n_vars * sizeof(*lp->cand_col_type));
    }
}


// TODO: Remove pop_size if there is no need to use.
void
allocate_population(MCproblem *mcp,  Population *pop, size_t pop_size)
{
    pop->size = pop_size;
    pop->indv = malloc(pop->size*sizeof(Individual));
    for (int i=0; i < pop->size; i++)
        allocate_individual(mcp, &(pop->indv[i]));
}

void
allocate_individual(MCproblem *mcp,  Individual *indv)
{
    indv->deletions = malloc(mcp->n_vars * sizeof(indv->deletions));
    if (mcp->use_modules)
        indv->modules = malloc(mcp->n_models * mcp->n_vars * sizeof(indv->modules));
    indv->objectives = malloc(mcp->n_models * sizeof(indv->objectives));
    indv->penalty_objectives = malloc(mcp->n_models * sizeof(indv->penalty_objectives));
}


void
free_population(MCproblem *mcp, Population *pop)
{
    for (int i=0; i < pop->size; i++)
        free_individual(mcp, &(pop->indv[i]));
    free(pop->indv);
}


void
free_individual(MCproblem *mcp, Individual *indv)
{
    free(indv->deletions);
    if (mcp->use_modules)
        free(indv->modules);
    free(indv->objectives);
    free(indv->penalty_objectives);
}


/* Sets individual variables randomly while  meeting constraints  */
void
set_random_individual(MCproblem *mcp,  Individual *indv)
{
    int i,j,k;
    int *deleted_rxns;
    SAFE_ALLOC(deleted_rxns = malloc(mcp->alpha * sizeof *deleted_rxns))

    /* init deletions */
    for (j = 0; j < mcp->n_vars; j++)
        indv->deletions[j] = !DELETED_RXN;
    for (i = 0; i < mcp->alpha; i++) {
        deleted_rxns[i] = (int)pcg32_boundedrand(mcp->n_vars);
        indv->deletions[deleted_rxns[i]] = DELETED_RXN;
    }
    /* init modules. Only one module reaction is inserted regardless of beta, this heuristic leads to better individuals */
    if (mcp->use_modules) {
        for (k = 0; k < mcp->n_models; k++) {
            for (j = 0; j < mcp->n_vars; j++)
                    indv->modules[k*mcp->n_vars + j] = !MODULE_RXN;
            indv->modules[k*mcp->n_vars + deleted_rxns[(int)pcg32_boundedrand(mcp->alpha)]] = MODULE_RXN;
        }
     }
    calculate_objectives(mcp, indv);
    free(deleted_rxns);
}

/* Individual without any genetic manipulations */
void
set_blank_individual(MCproblem *mcp,  Individual *indv)
{
    int j,k;
    /* init deletions */
    for (j = 0; j < mcp->n_vars; j++)
        indv->deletions[j] = !DELETED_RXN;
    /* init modules */
    if (mcp->use_modules) {
        for (k = 0; k < mcp->n_models; k++)
            for (j = 0; j < mcp->n_vars; j++)
                indv->modules[k*mcp->n_vars + j] = !MODULE_RXN;
     }
    for (k = 0; k < mcp->n_models; k++) {
        indv->objectives[k] = UNKNOWN_OBJ;
        indv->penalty_objectives[k] = UNKNOWN_OBJ;
    }
}


void
set_random_population(MCproblem *mcp, Population *pop)
{
    for (int i=0; i < pop->size; i++)
        set_random_individual(mcp, &(pop->indv[i]));
}

void
set_blank_population(MCproblem *mcp, Population *pop)
{
    for (int i=0; i < pop->size; i++)
        set_blank_individual(mcp, &(pop->indv[i]));
}
