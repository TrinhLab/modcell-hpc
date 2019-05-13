/* Core MOEA (NSGA-II) method, relays heavily on the methods defined in functions.c */

#include <stdlib.h>
#include "modcell.h"

void run_moea(MCproblem *mcp, Population *parent_population);
void selection_and_variation(MCproblem *mcp, Population *core_population, Population *offspring_population);
void evaluate_population(MCproblem *mcp, Population *population);
void environmental_selection(MCproblem *mcp, Population *parent_population, Population *offspring_population, Population *combined_population);
Individual * tournament_k2(MCproblem *mcp, Individual *indv1, Individual *indv2);

void assign_rank_and_crowding(MCproblem *mcp, Population *population);

/* Main loop of the MOEA
Notes:
- The combined population is allocated and copied, this is likely avoidable by representing it through pointers.
*/
void
run_moea(MCproblem *mcp, Population *parent_population)
{
    int n_generations = 0;
    double run_time = 0;
    Population *offspring_population = malloc(sizeof(Population));
    Population *combined_population = malloc(sizeof(Population));
    allocate_population(mcp, offspring_population, mcp->population_size);
    allocate_population(mcp, combined_population, 2*mcp->population_size);

    evaluate_population(mcp, parent_population);
    assign_rank_and_crowding(mcp, parent_population);
    n_generations++;

    while( (run_time < mcp->max_run_time) & (n_generations < mcp->n_generations) )
    {
        /* Core procedure */
        selection_and_variation(mcp, parent_population, offspring_population);
        evaluate_population(mcp, offspring_population);
        environmental_selection(mcp, parent_population, offspring_population, combined_population);

        /* Migration */

       /* Book keeping */
        n_generations++;
    }

    free_population(mcp, offspring_population);
    free_population(mcp, combined_population);
}


/* Creates offspring population by tournament selection, crossover, and mutation.
Notes:
- TODO: Candidate parents for tournament selection are selected purely at randomly. It might be valuable to consider a scheme where such candidates cannot repeat themselves, which might lead to better diversity.
- The new individuals will have some un-initialized fields.
*/
void
selection_and_variation(MCproblem *mcp, Population *parent_population, Population *offspring_population)
{
    int i;
    Individual *parent1, *parent2;

    /*Tournament selection and crossover*/
    for (i=0; i < mcp->population_size; i+=2){
        parent1 = tournament_k2(mcp, &(parent_population->indv[(int)pcg32_boundedrand(mcp->population_size)]), &(parent_population->indv[(int)pcg32_boundedrand(mcp->population_size)]));
        parent2 = tournament_k2(mcp, &(parent_population->indv[(int)pcg32_boundedrand(mcp->population_size)]), &(parent_population->indv[(int)pcg32_boundedrand(mcp->population_size)]));
        crossover(mcp, parent1, parent2, &(offspring_population->indv[i]), &(offspring_population->indv[i+1]));
    }

    /* Mutation and constraint enforcement */
    for (i=0; i < mcp->population_size; i++) {
        mutation(mcp, &(offspring_population->indv[i]));
        enforce_module_constraints(mcp, &(offspring_population->indv[i]));
    }
}


/* Tournament selection between two individuals for NSGA-II */
Individual *
tournament_k2(MCproblem *mcp, Individual *indv1, Individual *indv2)
{
    int f = find_domination(mcp, indv1, indv2);
    if (f == A_DOMINATES_B)  return (indv1);
    if (f == B_DOMINATES_A) return(indv2);
    if (indv1->crowding_distance > indv2->crowding_distance) return(indv1);
    if (indv2->crowding_distance > indv1->crowding_distance) return(indv2);
    return (pcg32_boundedrand(2) ? indv1 : indv2);
}


void
evaluate_population(MCproblem *mcp, Population *population)
{
    for (int i=0; i < mcp->population_size; i++)
        calculate_objectives(mcp, &(population->indv[i]));
}


void
assign_rank_and_crowding(MCproblem *mcp, Population *population)
{
    //PLACEHOLDA, instead, assign rank and crowding_distance
    Individual *indv;
    for (int i=0; i < mcp->population_size; i++){
        indv = &(population->indv[i]);
        indv->rank = 1;
        indv->crowding_distance = INF;
    }
}


/* Selects most fit individuals from both parents and offspring populations to create a new parent_population

Notes:
- Pass mixed population or pass separate populations? The second would not allow to call the function before an offspring has been created. This also relates to avoiding duplication in rank and crowsing calculations done in evaluate_population
*/
void
environmental_selection(MCproblem *mcp, Population *parent_population, Population *offspring_population, Population *combined_population)
{

    //PLACEHOLDA, fast non-dominated sorting will assign rank, crowding dist will be done separately.
    assign_rank_and_crowding(mcp, parent_population);

    // combine_populations(parent_population, offspring_population, combined_population)

    /* Non-dominated sorting */

    /* Assign crowing distance */

    /* Select the solutions in the last front based on their crowding distances */

    /* Modifiy core population */
}
