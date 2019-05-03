#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include "mcutils.h"
#include "modcell.h"

MCproblem readproblem(const char *problem_dir);
bool is_not_candidate(struct file ncandfile, const char *rxnid);
void load_parameters(MCproblem *mcp, const char *filepath);

MCproblem
read_problem(const char *problem_dir){
    MCproblem mcp = {NULL};
    DIR *d; //TODO: Move directory reading as an auxiliary function in mcutils.h
    struct dirent *dir;

    d = opendir(problem_dir);
    if (!d) {
        fprintf (stderr, "error: opening problem directory: '%s'.", problem_dir);
	exit(-1);
    }

    /* Obtain basic info */
    mcp.n_models = 0;
    while ((dir = readdir(d)) != NULL) {
        /* Determine number of models */
        if(strcmp( getExt (dir->d_name), ".mps") == 0){
    	    mcp.n_models++;
        }
        if(strcmp( dir->d_name, "cand") == 0){
            /* Determine deletion candidates */
    	    char full_path[256];
    	    strcpy(full_path, problem_dir);
    	    strcat(full_path, dir->d_name);
    	    struct file myfile = readFile(full_path);
            mcp.individual2id = myfile.array;
            mcp.n_vars = myfile.nlines;
    	    //freeFile(myfile); Cannot free since it contains information in mcp! See strdup()
    	}
    }
    closedir(d);

    /* Load problems */
    glp_smcp param;
    glp_init_smcp(&param);
    param.msg_lev = GLP_MSG_ERR;

    mcp.Ps = malloc(mcp.n_models * sizeof(*mcp.Ps));

    d = opendir(problem_dir);
    int i=0;
    while ((dir = readdir(d)) != NULL) {
        if(strcmp( getExt (dir->d_name), ".mps") == 0){
    	    char full_path[256];
    	    strcpy(full_path, problem_dir);
    	    strcat(full_path, dir->d_name);
            mcp.Ps[i] = glp_create_prob();
            glp_read_mps(mcp.Ps[i], GLP_MPS_FILE, NULL, full_path); // TODO: Silence this or redirect to log
            /* Calculate basis and solve so subsequent computations are warm-started */
            glp_adv_basis(mcp.Ps[i], 0);
            glp_simplex(mcp.Ps[i], &param);
            i++;
        }
    }
    closedir(d);

    /* Create maps between indvidual and bounds */
    mcp.individual2glp = malloc(mcp.n_models * sizeof(*mcp.individual2glp));
    for (int k=0; k < mcp.n_models; k++){
        mcp.individual2glp[k] = malloc(mcp.n_vars * sizeof(**mcp.individual2glp));

        char model_path[256];
    	strcpy(model_path, problem_dir);
    	strcat(strcat(model_path, glp_get_prob_name(mcp.Ps[k])), ".ncand");
    	struct file ncandfile = readFile(model_path);

      	glp_create_index(mcp.Ps[k]);
        for(int j=0; j < mcp.n_vars; j++){
            if(is_not_candidate(ncandfile, mcp.individual2id[j]))
                mcp.individual2glp[k][j] = NOT_CANDIDATE;
            else
                mcp.individual2glp[k][j] = glp_find_col(mcp.Ps[k], mcp.individual2id[j]); // If glp fails to find the column it will error
        }
      	glp_delete_index(mcp.Ps[k]);
    	freeFile(ncandfile);
    }

return mcp;
}

bool
is_not_candidate(struct file ncandfile, const char *rxnid){
    for (int i=0; i < ncandfile.nlines; i++)
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

    /* open file for reading   */
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

