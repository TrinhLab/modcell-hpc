#include <stdbool.h>
#include <glpk.h>
#include "pcg_basic.h"

#define NOT_CANDIDATE -1

typedef struct Individual{
 	/* ModCell */
 	bool *deletions; 		/* [n_variables]*/
     	bool **modules; 		/* [n_models][n_variables] */
	/* MOEA */
	double *objectives; 		/* [n_models] */
	double *penalty_objectives; 	/* [n_models] */
	int rank; // Should it be moved to pop?
	double crowding_distance;
} Individual;

typedef struct Population{
	Individual *indv;
	size_t size;
} Population;

typedef struct MCproblem{
	char *objective_type; // maybe use enum for this?
	unsigned int alpha;
	unsigned int beta;
	unsigned int n_models;
	char **model_names;

      	glp_prob **Ps;//[nmodels];
	int **individual2glp; //[nmodels][nvars] Contains model index that individual maps to or NOT_CANDIDATE if module is fixed.
	char **individual2id; //[nvars] Maps individual indices to reaction ID.
	// hash archive; // maps individuals deletions and module variables into objective values

	/* MOEA */
	unsigned int n_cores;
    	unsigned int seed;
    	int n_vars;
    	int population_size;
    	double pcross_bin;
    	double pmut_bin;
    	int n_generations;
    	int nbinmut;
    	int nbincross;
} MCproblem;


/* init.c */
void allocate_population(MCproblem *mcp, Population *indv, size_t size);
void free_population(MCproblem *mcp, Population *pop);
void allocate_individual(MCproblem *mcp, Individual *indv);
void free_individual(MCproblem *mcp, Individual *indv);
void set_random_population(MCproblem *mcp, Population *pop);
void set_blank_individual(MCproblem *mcp,  Individual *indv);
void set_random_individual(MCproblem *mcp,  Individual *indv);

/* functions.c */
double calculate_objectives(MCproblem *mcp, Individual *indv);

/* moea.c */
Population run_moea(MCproblem *mcp, Population *initial_population);
