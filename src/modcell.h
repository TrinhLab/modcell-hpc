#include <stdbool.h>
#include <glpk.h>
#include "pcg_basic.h"

#define NOT_CANDIDATE -1

//extern pcg32_random_t rng;

/* Main structures */

typedef struct {
	// modcell
	bool *deletions;//[nbits];
    	bool **modules;//[nmodels][nbits];
      	glp_prob **Ps;//[nmodels];
	// moea
	int rank;
	double crowding_distance;
	double *objectives;//[nmodels];
	double *penalty_objectives;//[nmodels];

} Individual;

//typedef struct
//{
//	individual *ind;
//} population;


// General problem information
typedef struct MCproblem
{
	char *objective_type; // maybe use enum for this?
	unsigned int alpha;
	unsigned int beta;
	unsigned int n_models;
      	glp_prob **Ps;//[nmodels];
	// refbounds *rbounds; This can be added later if glp query has too much overhead
	unsigned int **individual2glp; //[nmodels][nvars] Contains model index that individual maps to or NOT_CANDIDATE if module is fixed.
	char **individual2id; //[nvars] Maps individual indices to reaction ID.
	// hash archive; // maps individuals deletions and module variables into objective values

	//MOEA
	unsigned int n_cores;
    	unsigned int seed;
    	int n_vars;
//    int nobj; // called nmodels
    	int population_size;
    	double pcross_bin;
    	double pmut_bin;
    //double eta_c;
    //double eta_m;
    	int n_generations;
    int nbinmut;
    int nbincross;
    //int *nbits;
    //int bitlength;
    //int choice;
    //int obj1;
    //int obj2;
    //int obj3;
    //int angle1;
    //int angle2;
} MCproblem;


/* allocation.c */
Individual initialize_individual(MCproblem *mcp);
void free_individual(MCproblem *mcp, Individual *indv);
Individual **initialize_population(MCproblem *mcp);
void free_population(MCproblem *mcp, Individual *indv);

