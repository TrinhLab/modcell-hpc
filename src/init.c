/* Routines to allocate/free  and initialize objects*/

#include <stdlib.h>
#include "modcell.h"

void allocate_population(MCproblem *mcp,  Population *pop, size_t pop_size);
void allocate_individual(MCproblem *mcp,  Individual *indv);
void free_population(MCproblem *mcp, Population *pop);
void free_individual(MCproblem *mcp, Individual *indv);

void set_random_population(MCproblem *mcp, Population *pop);
void set_blank_individual(MCproblem *mcp,  Individual *indv);
void set_random_individual(MCproblem *mcp,  Individual *indv);

void
allocate_population(MCproblem *mcp,  Population *pop, size_t pop_size){
    pop->size = pop_size;
    pop->indv = (Individual *)malloc(pop->size*sizeof(Individual));
    for (int i=0; i<pop->size; i++)
        allocate_individual(mcp, &(pop->indv[i]));
}

//TODO: Make malloc pointer casting consistent
void
allocate_individual(MCproblem *mcp,  Individual *indv){
    indv->deletions = (bool *)malloc(mcp->n_vars * sizeof(bool));
    if (mcp->beta > 0){
        indv->modules = malloc(mcp->n_models * sizeof( *(indv->modules) ));
        for (int k = 0; k < mcp->n_models; k++)
            indv->modules[k] = malloc(mcp->n_vars * sizeof( **(indv->modules) ));
    }

        indv->objectives = (double *)malloc(mcp->n_models * sizeof(double));
        indv->penalty_objectives = (double *)malloc(mcp->n_models * sizeof(double));
}


void
free_population(MCproblem *mcp, Population *pop){
    for (int i=0; i<pop->size; i++)
        free_individual(mcp, &(pop->indv[i]));
    free(pop->indv);
}


void
free_individual(MCproblem *mcp, Individual *indv){
    free(indv->deletions);
    if (mcp->beta > 0){
        for (int k = 0; k < mcp->n_models; k++)
            free(indv->modules[k]);
        free(indv->modules);
    }
}


/* Sets individual variables randomly while  meeting constraints        */
void
set_random_individual(MCproblem *mcp,  Individual *indv){
    int i,j,k;
    int deleted_rxns[mcp->alpha]; // malloc?

    /* init deletions */
    for (j = 0; j < mcp->n_vars; j++)
        indv->deletions[j] = 1;
    for (i = 0; i < mcp->alpha; i++){
        deleted_rxns[i] = pcg32_boundedrand(mcp->n_vars);
        indv->deletions[deleted_rxns[i]] = 0;
    }
    /* init modules. Only one module reaction is inserted regardless of beta, this heuristic leads to better individuals*/
    if (mcp->beta > 0){
        for (int k = 0; k < mcp->n_models; k++){
            for (j = 0; j < mcp->n_vars; j++)
                    indv->modules[k][j] = 0;
            indv->modules[k][deleted_rxns[(int)pcg32_boundedrand(mcp->alpha)]] = 1;
        }
     }

    calculate_objectives(mcp, indv);
}

/* Individual without any genetic manipulations */
void
set_blank_individual(MCproblem *mcp,  Individual *indv){
    int i,j,k;
    /* init deletions */
    for (j = 0; j < mcp->n_vars; j++)
        indv->deletions[j] = 1;
    /* init modules */
    if (mcp->beta > 0){
        for (int k = 0; k < mcp->n_models; k++){
            for (j = 0; j < mcp->n_vars; j++)
                indv->modules[k][j] = 0;
        }
     }
    for (int k = 0; k < mcp->n_models; k++)
        indv->objectives[k] = -1; /* Use this to indivate objectives are not calculated */
}


void
set_random_population(MCproblem *mcp, Population *pop){
    for (int i=0; i < pop->size; i++)
        set_random_individual(mcp, &(pop->indv[i]));
}
