#include <stdbool.h>
#include <glpk.h>
#include "pcg_basic.h"

/* Definitions */
#define NOT_CANDIDATE -1

/* Globals */
glp_smcp param;

/* Structures */
typedef struct Individual{
 	/* ModCell */
 	bool *deletions; 		/* [n_variables]*/
     	bool **modules; 		/* [n_models][n_variables] */
	/* MOEA */ //IMPORTANT: Would these fields interfere or slow down the hash table look up? Maybe a hashkey struct should be created independently
	double *objectives; 		/* [n_models] */
	double *penalty_objectives; 	/* [n_models] */
	int rank;
	double crowding_distance;
} Individual;

typedef struct Population{
	Individual *indv;
	size_t size;
} Population;

typedef struct{
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

typedef struct MCproblem{
	char *objective_type; // maybe use enum for this or otherwise assert that the input is valid
	unsigned int alpha;
	unsigned int beta;
	unsigned int n_models;
	char **model_names;
	LPproblem *lps; /* [n_models] Contains everything needed to calculate an individuals fitness function */
	char **individual2id; /* [nvars] Maps individual indices to reaction ID. */
	/* MOEA */
	size_t n_cores;
    	size_t n_vars;
    	unsigned int population_size; //FIXME: setting this to size_t casuses seg fault
    	unsigned int seed;
    	unsigned int n_generations;
} MCproblem;


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

/* moea.c */
Population run_moea(MCproblem *mcp, Population *initial_population);
