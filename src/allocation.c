/* Routines to allocate and free objects*/

#include <stdlib.h>
#include "modcell.h"

void
allocate_population(MCproblem *mcp,  Population *pop, size_t pop_size){
    pop->size = pop_size;
    pop->indv = (Individual *)malloc(pop->size*sizeof(Individual));
    for (int i=0; i<pop->size; i++)
        allocate_individual(mcp, &(pop->indv[i]));
}

//TODO: Make malloc pointer casting consistent
void
allocate_individual_designvars(MCproblem *mcp,  Individual *indv){
    indv->deletions = (bool *)malloc(mcp->n_vars * sizeof(bool));
    if (mcp->beta > 0){
        indv->modules = malloc(mcp->n_models * sizeof( *(indv->modules) ));
        for (int k = 0; k < mcp->n_models; k++)
            indv->modules[k] = malloc(mcp->n_vars * sizeof( **(indv->modules) ));
    }
}

void
allocate_individual(MCproblem *mcp,  Individual *indv){
    allocate_individual_designvars(mcp, indv);
    /* New LP problem copies */
    indv->Ps = malloc(mcp->n_models * sizeof(*indv->Ps));
    for (int k = 0; k < mcp->n_models; k++){
        indv->Ps[k] = glp_create_prob();
        glp_copy_prob(indv->Ps[k], mcp->Ps[k], GLP_OFF); /* Note that trying to access any names will cause segfault due to GLP_OFF. */
    }
}

/* Allocate individual with external LP problems from another individual*/
void
partial_allocate_individual(MCproblem *mcp,  Individual *indv, Individual *old_indv){
    allocate_individual_designvars(mcp, indv);
    indv->Ps = old_indv->Ps; // TODO: test
}

void
free_population(MCproblem *mcp, Population *pop){
    for (int i=0; i<pop->size; i++)
        free_individual(mcp, &(pop->indv[i]));
    free(pop->indv);
}


/* Free individual without removing LP problems*/
void
partial_free_individual(MCproblem *mcp, Individual *indv){
    free(indv->deletions);
    if (mcp->beta > 0){
        for (int k = 0; k < mcp->n_models; k++)
            free(indv->modules[k]);
        free(indv->modules);
    }
}

void
free_individual(MCproblem *mcp, Individual *indv){
    partial_free_individual(mcp, indv);
    /* Free LPs */
    for (int k = 0; k < mcp->n_models; k++)
        glp_delete_prob(indv->Ps[k]);
    free(indv->Ps);
}




