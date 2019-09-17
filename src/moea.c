/* Core MOEA (NSGA-II) method, relays heavily on the methods defined in functions.c */

#include <stdlib.h>
#include <time.h>
#include "utlist.h"
#include "modcell.h"

void run_moea(MCproblem *mcp, Population *parent_population);
void selection_and_variation(MCproblem *mcp, Population *core_population, Population *offspring_population);
void evaluate_population(MCproblem *mcp, Population *population);
void environmental_selection(MCproblem *mcp, Population *parent_population, Population *offspring_population, Population *combined_population);
Individual * tournament_k2(MCproblem *mcp, Individual *indv1, Individual *indv2);
void add_individuals(MCproblem *mcp, Population *combined_pop, Population *parent_pop, item *head_fi, unsigned int fi_size, unsigned int *individuals_added);
void set_inf_crowding(MCproblem *mcp, Population *population);
void assign_crowding_distance(MCproblem *mcp, Population *pop, item *head_fi, unsigned int fi_size);

extern int mpi_pe, mpi_comm_size;

void migrate(MCproblem *mcp, Population *parent_population, Population *send_population, Population *receive_population);

/* Macros */
#define FREE_LIST(list_head) \
    DL_FOREACH_SAFE(list_head,elt,tmp) { \
      DL_DELETE(list_head,elt); \
      free(elt); \
    };

/*Function definitions */

/* Main loop of the MOEA
 * Notes:
 *      - The combined population is allocated and copied, this is likely avoidable by representing it through pointers.
 *      - When a PE finishes all others should terminate as well. This can be implemented with an MPI_Allreduce where one process sends 1. However, given the current synchronous nature of migration, syncing each generation would not add much of a performance overhad, and it is a much simple way of dealing with this problem.
 */
void
run_moea(MCproblem *mcp, Population *parent_population)
{
    unsigned int n_generations = 0;
    double run_time = 0;
    Population *offspring_population = malloc(sizeof(Population));
    Population *combined_population = malloc(sizeof(Population));
    allocate_population(mcp, offspring_population, mcp->population_size);
    allocate_population(mcp, combined_population, 2*mcp->population_size);

    Population *send_population = malloc(sizeof(Population));
    Population *receive_population = malloc(sizeof(Population));
    allocate_population(mcp, send_population, mcp->migration_size);
    allocate_population(mcp, receive_population, mcp->migration_size);

    /* Avoid uninitialized individuals  */
    set_blank_population(mcp, offspring_population);
    set_blank_population(mcp, combined_population);
    set_blank_population(mcp, send_population);
    set_blank_population(mcp, receive_population);

    clock_t begin = clock();
    evaluate_population(mcp, parent_population);
    n_generations++;

    set_inf_crowding(mcp, parent_population);
    set_inf_crowding(mcp, offspring_population);

    int done=0;
    while(1) {

        MPI_Barrier(MPI_COMM_WORLD); /* sync iterations. See notes in above function name for more details */
        if (run_time > mcp->max_run_time) { /* Simply re-check run time limit. E.g., one process has been waiting but was not done in prev iteration, but a slower one arrived here with done=1 due to time limit. Re-checking time limit should prevent issues. But mroe generally, when a process is done it should just send a message ot others to stop (specially to not attemt any further migrations)*/
            done = 1;
            if (mcp->verbose) printf("PE: %i\t Run time limit reached \t Time:%.1fs\n", mpi_pe, run_time);
        }
        if (done) break;

        /* Core procedure */
        selection_and_variation(mcp, parent_population, offspring_population);
        evaluate_population(mcp, offspring_population);
        environmental_selection(mcp, parent_population, offspring_population, combined_population);

        /* Migration */
        if ( (mpi_comm_size > 1) && (n_generations % mcp->migration_interval == 0) ) {
            if (mcp->verbose) printf("PE: %i Begin migration: %.0fs ...\n", mpi_pe, (double)(clock() - begin) / CLOCKS_PER_SEC);

            migrate(mcp, parent_population, send_population, receive_population);

            if (mcp->verbose) printf("...PE: %i end migration: %.0fs ...\n", mpi_pe, (double)(clock() - begin) / CLOCKS_PER_SEC);
        }

       /* Local book keeping */
        n_generations++;
        run_time = (double)(clock() - begin) / CLOCKS_PER_SEC;

        if (mcp->verbose && ( (n_generations-1) % PRINT_INTERVAL == 0))
            printf("PE: %i\t Generation:%i\t Time:%.1fs\n", mpi_pe, n_generations-1, run_time);

        if (run_time > mcp->max_run_time) {
            done = 1;
            if (mcp->verbose) printf("PE: %i\t Run time limit reached \t Time:%.1fs\n", mpi_pe, run_time);
        }
        if (n_generations > mcp->n_generations) {
            done = 1;
            if (mcp->verbose) printf("PE: %i\t Generation limit reached \t Time:%.1fs\n", mpi_pe, run_time);
        }
    }

    free_population(mcp, offspring_population);
    free_population(mcp, combined_population);
    free_population(mcp, send_population);
    free_population(mcp, receive_population);
}



/* Sends and receives individuals from other islands (PEs)
 * Notes:
 *      - Migration policies other than random depend on how parent_population is sorted.
 *      - MPI is designed to send arrays of structures made of simple types, but not structure of arrays. So the fields of the individual struct are sent independently, which of course requires more messages.
 * Asynchronous notes:
 *      - Individuals could  be replaced directly in the parent population, but copying to send and receive populations allows asynchronous message passing. It can also avoid data race conditions.
 *      - Currently non-blocking communication is likely not beneficial. To increase performance with non-blocking message passing, use migration flag in main loop that can be checked each generation, migrate sets it to 0 if it succeeds to complete both send and receive.
 *      - Using async calls also allows for a timeout mechanism or some other way to cancel the message passing. This might be necessary since a PE might be in the middle of sending migration messages while another one is finished? The barrier (coupled with other PE termination message) should prevent that.
 *      - FIXME: Current blocking calls could end up in deadlock since messages "cross" each other. So this should use non-blocking messages.
 */

void
migrate(MCproblem *mcp, Population *parent_population, Population *send_population, Population *receive_population)
{
    int target_pe;
    int i, j, k, tag = 10;
    Individual *send_indv, *recv_indv;
    MPI_Status status;
    int *send_idx = malloc(mcp->migration_size * sizeof(int));
    int *receive_idx = malloc(mcp->migration_size * sizeof(int));

    /* Message passing topology */ // If it is not dynamic it can be determined outside of this method
    if (mcp->migration_topology == MIGRATION_TOPOLOGY_RING) {
        if (mpi_pe == mpi_comm_size - 1)
            target_pe = 0;
        else
            target_pe = mpi_pe + 1;
    }
    else if (mcp->migration_topology == MIGRATION_TOPOLOGY_RANDOM) {
        do {
            target_pe = (int)pcg32_boundedrand(mpi_comm_size);
        } while (target_pe != mpi_pe);
    } else { fprintf (stderr, "error: Invalid migration topology option"); exit(-1); }

    /* Determine individuals to send and receive */
    if (mcp->migration_policy == MIGRATION_POLICY_REPLACE_SENT) {
        for (i=0; i < mcp->migration_size; i++) {
            send_idx[i] = i;
            receive_idx[i] = i;
        }
    }
    else if (mcp->migration_policy == MIGRATION_POLICY_REPLACE_BOTTOM) {
        for (i=0; i < mcp->migration_size; i++) {
            send_idx[i] = i;
            receive_idx[i] = mcp->population_size - 1 - i;
        }
    }
    else if (mcp->migration_policy == MIGRATION_POLICY_RANDOM) {
        for (i=0; i < mcp->migration_size; i++) {
            send_idx[i] = i; (int)pcg32_boundedrand(mcp->population_size);
            receive_idx[i] = send_idx[i];
        }
    } else { fprintf (stderr, "error: Invalid migration policy option"); exit(-1); }


    /* Message exchange */
    for (i=0; i < mcp->migration_size; i++) {

        copy_individual(mcp,  &(parent_population->indv[send_idx[i]]), &(send_population->indv[i]));
        send_indv = &(send_population->indv[i]);
        recv_indv = &(receive_population->indv[i]);

        /* Deletions */
        MPI_Send(send_indv->deletions, mcp->n_vars, MPI_INT, target_pe, tag, MPI_COMM_WORLD);
        MPI_Recv(recv_indv->deletions, mcp->n_vars, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        /* Module reactions */
        if (mcp->beta > 0) {
            for (k=0; k < mcp->n_models; k++) {
                // Sending the array leads to segfaults in glpk...
                //MPI_Send(&(send_indv->modules[k]), mcp->n_vars, MPI_INT, target_pe, tag, MPI_COMM_WORLD);
                //MPI_Recv(&(recv_indv->modules[k]), mcp->n_vars, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                for (j=0; j < mcp->n_vars; j++) {
                    MPI_Send(&(send_indv->modules[k][j]), 1, MPI_INT, target_pe, tag, MPI_COMM_WORLD);
                    MPI_Recv(&(recv_indv->modules[k][j]), 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                }
            }
        }
        /* Objectives */
        MPI_Send(send_indv->objectives, mcp->n_models, MPI_DOUBLE, target_pe, tag, MPI_COMM_WORLD);
        MPI_Recv(recv_indv->objectives, mcp->n_models, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        /* Penalty_Objectives */
        MPI_Send(send_indv->penalty_objectives, mcp->n_models, MPI_DOUBLE, target_pe, tag, MPI_COMM_WORLD);
        MPI_Recv(recv_indv->penalty_objectives, mcp->n_models, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        /* Crowding distance */
        MPI_Send(&(send_indv->crowding_distance), 1, MPI_DOUBLE, target_pe, tag, MPI_COMM_WORLD);
        MPI_Recv(&(recv_indv->crowding_distance), 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    /* Copy received individuals into parent population */
    for (i=0; i < mcp->migration_size; i++)
        copy_individual(mcp,  &(receive_population->indv[i]), &(parent_population->indv[receive_idx[i]]));

    free(send_idx);
    free(receive_idx);

    MPI_Barrier(MPI_COMM_WORLD); // This is for rare cases and general safety. Specially for certain topologies.E.g., one process finished and writes its population while an others tries to send to it. Given the current blocking calls the performance overhead is likely irrelevant. If async is implemented this obviously needs to be modified.
}


/* Creates offspring population by tournament selection, crossover, and mutation.
 * Notes:
 *      - TODO: Candidate parents for tournament selection are selected purely at random. It might be valuable to consider a scheme where such candidates cannot repeat themselves, which might lead to better diversity.
 *      - The new individuals will have some uninitialized fields.
 */
void
selection_and_variation(MCproblem *mcp, Population *parent_population, Population *offspring_population)
{
    int i;
    Individual *parent1, *parent2;

    /*Tournament selection and crossover*/
    for (i=0; i < mcp->population_size; i+=2) {
        parent1 = tournament_k2(mcp, &(parent_population->indv[(int)pcg32_boundedrand(mcp->population_size)]), &(parent_population->indv[(int)pcg32_boundedrand(mcp->population_size)]));
        parent2 = tournament_k2(mcp, &(parent_population->indv[(int)pcg32_boundedrand(mcp->population_size)]), &(parent_population->indv[(int)pcg32_boundedrand(mcp->population_size)]));
        crossover(mcp, parent1, parent2, &(offspring_population->indv[i]), &(offspring_population->indv[i+1]));
    }

    /* Mutation */
    for (i=0; i < mcp->population_size; i++)
        mutation(mcp, &(offspring_population->indv[i]));

    /* Constraint enforcement */
    if (mcp->beta > 0)
        for (i=0; i < mcp->population_size; i++)
            enforce_module_constraints(mcp, &(offspring_population->indv[i]));
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

/* Assign objective values to each individual */
void
evaluate_population(MCproblem *mcp, Population *population)
{
    for (int i=0; i < mcp->population_size; i++)
        calculate_objectives(mcp, &(population->indv[i]));
}

/* Selects most fit individuals from both parents and offspring populations to create a new parent_population
 * Notes:
 *      - Do non-dominated sorting of each front, calculate distances (if needed, i.e., front size is greater than individuals left to fill pop), then stop when pop is filled
 */
void
environmental_selection(MCproblem *mcp, Population *parent_pop, Population *offspring_pop, Population *combined_pop)
{
    combine_populations(mcp, parent_pop, offspring_pop, combined_pop);

    int i, it1, it2, domflag, is_fi_empty;
    int *domination_count = calloc(combined_pop->size, sizeof(int)); /* keeps domination count (np, nq) of each individual */
    Individual *p, *q;
    item *head_sp = NULL, *head_fi = NULL, *head_f_next=NULL, *node, *tmp, *elt, *ref, *elt_p, *elt_q;
    unsigned int individuals_added = 0, fi_size = 0;

    for (it1=0; it1 < combined_pop->size; it1++) {
        for (it2=0; it2 < combined_pop->size; it2++) {
            p = &(combined_pop->indv[it1]);
            q = &(combined_pop->indv[it2]);

            domflag = find_domination(mcp, p, q);
            if (domflag == A_DOMINATES_B) {
                SAFE_ALLOC(node = (item *)malloc(sizeof *node))
                node->index = it2;
                DL_APPEND(head_sp, node); /* add q to Sp */
            } else if (domflag == B_DOMINATES_A) {
                domination_count[it1]++;
            }
        }
        if (domination_count[it1] == 0) {
            p->rank = 1;
            SAFE_ALLOC(node = (item *)malloc(sizeof *node))
            node->index = it1;
            DL_APPEND(head_fi, node); /* add p to F1 */
            fi_size++;
        }
    }

    add_individuals(mcp, combined_pop, parent_pop, head_fi, fi_size, &individuals_added);
    if(individuals_added == parent_pop->size)
        goto cleanup;

    i = 1; /* Fornt N. */
    is_fi_empty = 0;
    head_f_next = NULL;
    while (!is_fi_empty) {
        is_fi_empty = 1;
        DL_FOREACH(head_fi, elt_p) {
            DL_FOREACH_SAFE(head_sp, elt_q, tmp) {
                p = &(combined_pop->indv[elt_p->index]);
                q = &(combined_pop->indv[elt_q->index]);

                domination_count[elt_q->index]--;
                if (domination_count[elt_q->index] == 0) {
                    q->rank = i+1;
                    DL_DELETE(head_sp, elt_q);
                    DL_APPEND(head_f_next, elt_q);
                    is_fi_empty=0;
                }
            }
        }
        i +=1;

        /* Overwrite fi with f_next */
        FREE_LIST(head_fi)
        fi_size = 0;
        ref = head_fi;
        DL_FOREACH_SAFE(head_f_next, elt, tmp) {
            DL_DELETE(head_f_next, elt);
            DL_APPEND_ELEM(head_fi, ref, elt);
            ref = elt;
            fi_size++;
        }

        add_individuals(mcp, combined_pop, parent_pop, head_fi, fi_size, &individuals_added);
        if(individuals_added == parent_pop->size)
            goto cleanup;
    }

    /* Free memory */
cleanup:   free(domination_count);
    FREE_LIST(head_sp);
    FREE_LIST(head_fi);
    FREE_LIST(head_f_next);
}



/* UTLIST: "The comparison function must return an int that is negative, zero, or positive, which specifies whether the first item should sort before, equal to, or after the second item, respectively." (i.e., sort is in ascending order)
*/
Population *g_pop;
int g_k;
int
sortcmp(item *a, item *b)
{
    if (g_pop->indv[a->index].penalty_objectives[g_k] > g_pop->indv[b->index].penalty_objectives[g_k])
        return -1;
    if (g_pop->indv[a->index].penalty_objectives[g_k] < g_pop->indv[b->index].penalty_objectives[g_k])
        return 1;
    //else if (g_pop->indv[a->index].penalty_objectives[g_k] == g_pop->indv[b->index].penalty_objectives[g_k])
    return 0;
}

/* Notes:
 *   - A greater value of crowding distance is better, since the edge individuals to be preserved obtain a crowding distance of INF.
 */
int
sortcrowding(item *a, item *b)
{
    if (g_pop->indv[a->index].crowding_distance > g_pop->indv[b->index].crowding_distance)
        return -1;
    if (g_pop->indv[a->index].crowding_distance < g_pop->indv[b->index].crowding_distance)
        return  1;
    return 0;
}

/* Adds individuals and determines if crowding distance needs to be calculated
 *      - list_head is the list containing a  non-dominated front.
 *      - last_index is the index of the last individual in new_pop.
 * Notes:
 *      - Sorting a list pointer within a function causes errors when trying to free that list. The solution is to work with local copies (or to use a library other than UTLIST). That is also why this function has the crowding_distance calculation embedded instead of in a smaller function
 *      - Alternative metrics can be used instead of crowding distance, such as reference point distance.
 */

void
add_individuals(MCproblem *mcp, Population *combined_pop, Population *parent_pop, item *head_fi_og, unsigned int fi_size, unsigned int *individuals_added)
{
    assert(*individuals_added < parent_pop->size);

    item *head_fi = NULL, *elt, *node, *tmp;

    /* Create local copy of head_fi */
    DL_FOREACH(head_fi_og, elt) {
        SAFE_ALLOC(node = (item *)malloc(sizeof *node))
        node->index = elt->index;
        DL_APPEND(head_fi, node);
    }

    item *elt_p;

    /* If fi size is above what is needed to fill the pop calculate crowding distance and sort fi by it*/

    if (fi_size > (parent_pop->size - *individuals_added)) {
        item  *elt, *tail_fi;
        int it;
        double fm_max, fm_min;

        DL_FOREACH(head_fi, elt) {
            combined_pop->indv[elt->index].crowding_distance = 0;
        }

        g_pop = combined_pop; /* Passes population to sort funcitons*/
        for(int k=0; k < mcp->n_models; k++) {
            g_k = k; /* Passes objective index to sortcmp() */

            DL_SORT(head_fi, sortcmp);
            DL_FOREACH(head_fi, elt) { /* get tail of the list */
                tail_fi = elt;
            }

            combined_pop->indv[head_fi->index].crowding_distance = INF;
            combined_pop->indv[tail_fi->index].crowding_distance = INF;

            fm_max = combined_pop->indv[head_fi->index].penalty_objectives[k];
            fm_min = combined_pop->indv[tail_fi->index].penalty_objectives[k];

            if (fm_max != fm_min) { /* Avoid calculating c.d. with 0 division */
                it = 0;
                DL_FOREACH(head_fi, elt) {
                    if ( (it > 0) && (it < fi_size -1) ){
                        combined_pop->indv[elt->index].crowding_distance = combined_pop->indv[elt->index].crowding_distance
                            + ( combined_pop->indv[elt->next->index].penalty_objectives[k]
                            - combined_pop->indv[elt->prev->index].penalty_objectives[k] )
                            / (fm_max - fm_min);
                    }
                    it++;
                }
            }
        }

        /* Sort fi_size by crowding distance */
        DL_SORT(head_fi, sortcrowding);
    }

    DL_FOREACH(head_fi, elt_p) {
        copy_individual(mcp,  &(combined_pop->indv[elt_p->index]), &(parent_pop->indv[*individuals_added]));
        *(individuals_added) +=  1;
        if (*individuals_added == parent_pop->size)
            return;
    }

    FREE_LIST(head_fi); /* free local copy of head_fi */

}

/* Makes sure that crowding distance is assigned for tournament selection */
void
set_inf_crowding(MCproblem *mcp, Population *population)
{
    for (int i=0; i < mcp->population_size; i++){
        population->indv[i].crowding_distance = INF;
    }
}

