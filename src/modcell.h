#include <stdbool.h>
#include <stdio.h>
#include <glpk.h>
#include "pcg_basic.h"
#include "mpi.h"

/* Notation */
#define UNKNOWN_OBJ -1
#define NOT_CANDIDATE -1
#define DELETED_RXN 0
#define MODULE_RXN 1
#define A_DOMINATES_B 1
#define B_DOMINATES_A -1
#define NONDOMINATED 0

/* Parameter notation */
#define MIGRATION_POLICY_REPLACE_BOTTOM 0
#define MIGRATION_POLICY_REPLACE_SENT 1
#define MIGRATION_POLICY_RANDOM 2
#define MIGRATION_TOPOLOGY_RING 0
#define MIGRATION_TOPOLOGY_RANDOM 1

/* Definitions */
#define INF 1.0e14 		/* A value to simulate infinity */
#define MAX_MODULES 200 	/* A value above any practical beta expected, used for array allocation. FIXME: This should be done dynamically*/
#define LP_TIME_LIMIT_MILISEC 10000 /* Ensures GLPK does not get stuck trying to solve an LP */
#define LP_MSG_LEV GLP_MSG_OFF 	/* GLP output, options are: GLP_MSG_ERR  (will sometimes indicate that an LP could not be solved due to numerical issues), GLP_MSG_ALL (usefull for debuggin), or GLP_MSG_OFF (to avoid output)*/
#define OBJ_TOL 0.015 		/* Tolerance value to consider two objectives different. Currently only look at two decimal digits, the 0.005 in the last place is for rounding */

/* Parameters */
#define PRINT_INTERVAL 10 	/* Generations interval when info is printed */

/* Macros */
#define SAFE_ALLOC(expr) if( (expr) == NULL) { printf("Memory allocation failed, exiting...\n"); exit(-1);}

/* ifdef settings */
#define MIN_LOG 0 		/* Use it to work around GLPK un-silenceable output. However, turning this on messes up output buffering in MPI so only PE=0 prints in real time, while the rest print at the end */

/* Globals */
glp_smcp param;
int mpi_pe, mpi_comm_size;

/* Structures */

typedef struct {
 	/* ModCell */
 	bool *deletions; 		/* [n_variables]*/
     	bool *modules; 			/* [n_models*n_variables] */
	/* MOEA */
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
	/* modcell */
	char objective_type[256]; // TODO use enum for this or otherwise assert that the input is valid
	unsigned int alpha;
	unsigned int beta;
	unsigned int n_models;
	char **model_names; 	/* [nvars] */
	LPproblem *lps; 	/* [n_models] Contains everything needed to calculate an individuals fitness function */
	char **individual2id; 	/* [nvars] Maps individual indices to reaction ID. */

	/* MOEA */
    	size_t n_vars;
    	unsigned int population_size; // TODO: Use size_t consistently
    	unsigned int seed; /* Note: The real RNG seed is seed + MPI PE number */
    	unsigned int n_generations;
	double crossover_probability;
	double mutation_probability;
	double max_run_time; 	/* maximum run time in seconds */

	/* Parallelization  */
    	unsigned int migration_interval;
    	unsigned int migration_size;
    	unsigned int migration_policy;
    	unsigned int migration_topology;

	/* Other */
	int verbose;
    	int use_modules;  /* = hmcp.beta > 0 */
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
void set_random_individual(MCproblem *mcp,  Individual *indv);
void set_blank_population(MCproblem *mcp, Population *pop);
void set_blank_individual(MCproblem *mcp,  Individual *indv);

/* functions.c */
void calculate_objectives(MCproblem *mcp, Individual *indv);
void crossover(MCproblem *mcp, Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
void mutation(MCproblem *mcp, Individual *indv);
void enforce_module_constraints(MCproblem *mcp, Individual *indv);
int find_domination(MCproblem *mcp, Individual *indv_a, Individual *indv_b);
void copy_individual(MCproblem *mcp, Individual *indv_source, Individual *indv_dest);
void combine_populations(MCproblem *mcp, Population *pop1, Population *pop2, Population *combined_pop);
void calculate_objective(MCproblem *mcp, Individual *indv, int k, int *change_bound);

/* moea.c */
void run_moea(MCproblem *mcp, Population *initial_population);

/* module_minimizer.c */
void minimize_mr(MCproblem *mcp, Population *parent_population);
