#define NOT_CANDIDATE -1

#include <stdbool.h>
#include <glpk.h>

typedef struct
{
	// modcell
	bool *deletions;//[nbits];
    	bool *modules;//[nmodels][nbits];
      	glp_prob **P;//[nmodels];
	// moea
	int rank;
	double crowd_dist;
	double *obj;//[nmodels];
	double *penaltyobj;//[nmodels];

} individual;

typedef struct
{
	individual *ind;
} population;


// General problem information
typedef struct MCprob
{
	char *objective;
	int alpha;
	int beta;
	int nmodels;
      	glp_prob **Ps;//[nmodels];
	// refbounds *rbounds; This can be added later if glp query has too much overhead
	int **individual2glp; //[nmodels][nvars] Contains model index that individual maps to or NOT_CANDIDATE if module is fixed.
	char **individual2id; //[nvars] Maps individual indices to reaction ID.
	// hash archive; // maps individuals deletions and module variables into objective values

	//MOEA
	int ncores;
    double seed;
    int nvars;
//    int nobj; // called nmodels
    int popsize;
    double pcross_bin;
    double pmut_bin;
    //double eta_c;
    //double eta_m;
    int ngen;
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
} MCprob;

