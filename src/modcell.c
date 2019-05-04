#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include "mcutils.h"
#include "modcell.h"


// Move to header? Should modcell.h be renamed? e.g. to general.h?
MCproblem read_problem(const char *problem_dir);
bool is_not_candidate(Charlist ncandfile, const char *rxnid);
void load_parameters(MCproblem *mcp, const char *filepath);


MCproblem
read_problem(const char *problem_dir_path){
    MCproblem mcp = {NULL};
    int i,k,j;

    Charlist pd = read_dir(problem_dir_path);

    /* Obtain basic info */
    mcp.n_models = 0;
    for (i=0; i < pd.n; i++){
        /* Determine number of models */
        if(is_extension(pd.array[i], "mps")){
    	    mcp.n_models++;
        }
        /* Determine deletion candidates */
        if(is_extension(pd.array[i], "cand")){
            char cand_path[256];
            set_full_path(cand_path, problem_dir_path, pd.array[i]);
    	    Charlist f = read_file(cand_path);
            mcp.n_vars = f.n;
            mcp.individual2id = malloc(f.n * sizeof *mcp.individual2id);
            for (i=0; i < f.n; i++)
                mcp.individual2id[i] = strdup(f.array[i]);
    	    free_charlist(f);
    	}
    }

    /* Load problems */
    glp_smcp param;
    glp_init_smcp(&param);
    param.msg_lev = GLP_MSG_ERR;

    mcp.Ps = malloc(mcp.n_models * sizeof(*mcp.Ps));
    k=0;
    for (i=0; i < pd.n; i++){
        if(is_extension(pd.array[i], "mps")){
            char prob_path[256];
            set_full_path(prob_path, problem_dir_path, pd.array[i]);
            mcp.Ps[k] = glp_create_prob();
            glp_read_mps(mcp.Ps[k], GLP_MPS_FILE, NULL, prob_path); // TODO: Silence this or redirect to log
            /* Calculate basis and solve so subsequent computations are warm-started */
            glp_adv_basis(mcp.Ps[k], 0);
            glp_simplex(mcp.Ps[k], &param);
            k++;
        }
    }

    free_charlist(pd);

    /* Create maps between indvidual and bounds */
    mcp.individual2glp = malloc(mcp.n_models * sizeof(*mcp.individual2glp));
    for (k=0; k < mcp.n_models; k++){
        mcp.individual2glp[k] = malloc(mcp.n_vars * sizeof(**mcp.individual2glp));

        char model_path[256];
    	strcpy(model_path, problem_dir_path);
    	strcat(strcat(model_path, glp_get_prob_name(mcp.Ps[k])), ".ncand");
    	Charlist ncandfile = read_file(model_path);

      	glp_create_index(mcp.Ps[k]);
        for(j=0; j < mcp.n_vars; j++){
            if(is_not_candidate(ncandfile, mcp.individual2id[j]))
                mcp.individual2glp[k][j] = NOT_CANDIDATE;
            else
                mcp.individual2glp[k][j] = glp_find_col(mcp.Ps[k], mcp.individual2id[j]); /* If glp fails to find the column it will error */
        }
      	glp_delete_index(mcp.Ps[k]);
    	free_charlist(ncandfile);
    }

return mcp;
}

bool
is_not_candidate(Charlist ncandfile, const char *rxnid){
    for (int i=0; i < ncandfile.n; i++)
        if(!strcmp(ncandfile.array[i], rxnid))
                return true;
    return false;
}

#define READPARAM(id, type) \
    if (fscanf(fp, #id"=%"#type"\n", &mcp->id) !=1){ \
        fprintf (stderr, "error: parsing parameter'%s'\n", #id); \
	exit(-1); \
    } // print parameter here if interested

void
load_parameters(MCproblem *mcp, const char *filepath){
    FILE *fp = NULL;

    if (!(fp = fopen (filepath, "r"))) {
        fprintf (stderr, "error: file open failed '%s'\n", filepath);
	exit(-1);
    }
    READPARAM(objective_type,s);
    READPARAM(alpha,d);
    READPARAM(beta,d);
    READPARAM(population_size,d);
    READPARAM(n_generations,d);
    READPARAM(n_cores,d);
    READPARAM(seed,u);
    if (fp) fclose (fp);
}


//void
//writepopulation(MCproblem *mcp, Individual *outpop, const char *out_population_path){
//}
//
//void
//set_initial_population(MCproblem *mcp, Individual **pop, char *population_path){
//};
//
//void
//set_random_population(MCproblem *mcp, Individual **pop){
//}

// population runmoea(MCproblem *mcp) // define on separate file

int
main (int argc, char **argv) {

/* Parse  input */
if (argc < 3){
	fprintf (stderr, "error: insufficient input, usage: %s <problem directory> <parameters> <initial population (optional)>\n", argv[0]);
        return(1);
    }

MCproblem mcp = read_problem(argv[1]);
load_parameters(&mcp, argv[2]);
printf("--------------------------------------------------\n");

/* Seed global RNG */
pcg32_srandom(mcp.seed, 54u);

/* Intialize population */
// allocate
if (argc == 4){
    char pop_path[256];
    printf("Initial population path: %s\n", argv[3]);
    strcpy(pop_path, argv[3]);
    // set population to file read
}
else{
    printf("Initial population not specified\n");
    // set random population
}

/* Run */

return(0);
}

