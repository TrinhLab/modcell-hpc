#include <stdbool.h>
#include <glpk.h>
#include "pcg_basic.h"

/* Notation */
#define NOT_CANDIDATE -1
#define DELETED_RXN 0
#define MODULE_RXN 1
#define A_DOMINATES_B 1
#define B_DOMINATES_A -1
#define NONDOMINATED 0

/* Definitions */
#define INF 1.0e14 /* A value to simulate infinity */
#define MAX_MODULES 200 /* A value above any practical beta expected, used for array allocation. */

/* Macros */
#define SAFE_ALLOC(expr) if( (expr) == NULL) { printf("Memory allocation failed, exiting...\n"); exit(-1);}

/* Globals */
glp_smcp param;

/* Structures */
typedef struct {
 	/* ModCell */
 	bool *deletions; 		/* [n_variables]*/
     	bool **modules; 		/* [n_models][n_variables] */
	/* MOEA */ //IMPORTANT: Would these fields interfere (or slow down) the hash table look up? Maybe a hashkey struct should be created independently
	double *objectives; 		/* [n_models] */
	double *penalty_objectives; 	/* [n_models] */
	int rank; // This is currently unused.
	double crowding_distance;
} Individual;

typedef struct Population{
	Individual *indv;
	size_t size;
} Population;

typedef struct {
	// char *name; /* This is already in model_names but might be useful for bookeeping */
	glp_prob *P; 		/* GLPK LP problem */
	int *cand_col_idx; 	/* [nvars] Contains model index that individual maps to or NOT_CANDIDATE if module is fixed. */
	double *cand_og_lb; 	/* [n_cands] Maps indices of individuals to original lower bound values */
	double *cand_og_ub; 	/* [n_cands] Maps indices of individuals to original lower bound values */
	int *cand_col_type; 	/* [n_vars] GLPK column type */
	int prod_col_idx; 	/* Index of the product secretion reaction */
	int bio_col_idx; 	/* Index of the biomass formation reaction */
	double max_prod_growth; /* Maximum rate of product synthesis for growth state */
	double no_deletion_objective; /* Objective value when no deletions are present */
} LPproblem;

typedef struct {
	char objective_type[256]; // maybe use enum for this or otherwise assert that the input is valid
	unsigned int alpha;
	unsigned int beta;
	unsigned int n_models;
	char **model_names; 	/* [nvars] */
	LPproblem *lps; 	/* [n_models] Contains everything needed to calculate an individuals fitness function */
	char **individual2id; 	/* [nvars] Maps individual indices to reaction ID. */
	/* MOEA */
    	size_t n_vars;
	unsigned int n_cores;
    	unsigned int population_size; //FIXME: Use size_t consistently
    	unsigned int seed;
    	unsigned int n_generations;
	double crossover_probability;
	double mutation_probability;
	double max_run_time;
} MCproblem;

typedef struct item { /* list item */
     int index;
     struct item *prev, *next;
} item;

/* init.c */
void allocate_MCproblem(MCproblem *mcp, unsigned int n_models, size_t n_vars);
void allocate_population(MCproblem *mcp, Population *indv, size_t size);
void free_population(MCproblem *mcp, Population *pop);
void allocate_individual(MCproblem *mcp, Individual *indv);
void free_individual(MCproblem *mcp, Individual *indv);
void set_random_population(MCproblem *mcp, Population *pop);
void set_blank_individual(MCproblem *mcp,  Individual *indv);
void set_random_individual(MCproblem *mcp,  Individual *indv);

/* functions.c */
void calculate_objectives(MCproblem *mcp, Individual *indv);
void crossover(MCproblem *mcp, Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
void mutation(MCproblem *mcp, Individual *indv);
void enforce_module_constraints(MCproblem *mcp, Individual *indv);
int find_domination(MCproblem *mcp, Individual *indv_a, Individual *indv_b);
void copy_individual(MCproblem *mcp, Individual *indv_source, Individual *indv_dest);
void combine_populations(MCproblem *mcp, Population *pop1, Population *pop2, Population *combined_pop);

/* moea.c */
void run_moea(MCproblem *mcp, Population *initial_population);
