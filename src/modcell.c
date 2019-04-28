#include <stdio.h>
#include <string.h>
#include "mcutils.h"
#include "modcell.h"
#include <dirent.h>

MCprob readproblem(const char *problemdir);
bool isnotcandidate(struct file ncandfile, const char *rxnid);
void loadparams(MCprob *mcprob, const char *filepath);

MCprob
readproblem(const char *problemdir){
    MCprob mcprob = {NULL};
    DIR *d;
    struct dirent *dir;

    d = opendir(problemdir);
    if (!d) {
        fprintf (stderr, "error: opening problem directory: '%s'.", problemdir);
	exit(-1);
    }

    /* Obtain basic info */
    mcprob.nmodels = 0;
    while ((dir = readdir(d)) != NULL) {
        /* Determine number of models */
        if(strcmp( getExt (dir->d_name), ".mps") == 0){
    	    mcprob.nmodels++;
        }
        if(strcmp( dir->d_name, "cand") == 0){
            /* Determine deletion candidates */
    	    char full_path[256];
    	    strcpy(full_path, problemdir);
    	    strcat(full_path, dir->d_name);
    	    struct file myfile = readFile(full_path);
            mcprob.individual2id = myfile.array;
            mcprob.nvars = myfile.nlines;
    	    //freeFile(myfile); Cannot free since it contains information in mcprob!
    	}
    }
    closedir(d);

    /* Load problems */
    glp_smcp param;
    glp_init_smcp(&param);
    param.msg_lev = GLP_MSG_ERR;

    mcprob.Ps = malloc(mcprob.nmodels * sizeof(*mcprob.Ps)); // TODO: why does changing nmodels to lower value still work?

    d = opendir(problemdir);
    int i=0;
    while ((dir = readdir(d)) != NULL) {
        if(strcmp( getExt (dir->d_name), ".mps") == 0){
    	    char full_path[256];
    	    strcpy(full_path, problemdir);
    	    strcat(full_path, dir->d_name);
            mcprob.Ps[i] = glp_create_prob();
            glp_read_mps(mcprob.Ps[i], GLP_MPS_FILE, NULL, full_path); // TODO: Silence this or redirect to log
            /* Calculate basis and solve so subsequent computations are warm-started */
            glp_adv_basis(mcprob.Ps[i], 0);
            glp_simplex(mcprob.Ps[i], &param);
            i++;
        }
    }
    closedir(d);

    /* Create maps between indvidual and bounds */
    mcprob.individual2glp = malloc(mcprob.nmodels * sizeof(*mcprob.individual2glp));
    for (int k =0; k < mcprob.nmodels; k++){
        mcprob.individual2glp[k] = malloc(mcprob.nvars * sizeof(**mcprob.individual2glp));

        char model_path[256];
    	strcpy(model_path, problemdir);
    	strcat(strcat(model_path, glp_get_prob_name(mcprob.Ps[k])), ".ncand");
    	struct file ncandfile = readFile(model_path);

      	glp_create_index(mcprob.Ps[k]);
        for(int j=0; j < mcprob.nvars; j++){
            if(isnotcandidate(ncandfile, mcprob.individual2id[j]))
                mcprob.individual2glp[k][j] = NOT_CANDIDATE;
            else
                mcprob.individual2glp[k][j] = glp_find_col(mcprob.Ps[k], mcprob.individual2id[j]); // If glp fails to find the column it will error
        }
      	glp_delete_index(mcprob.Ps[k]);
    	freeFile(ncandfile);
    }

return mcprob;
}

bool
isnotcandidate(struct file ncandfile, const char *rxnid){
    for (int i=0; i < ncandfile.nlines; i++)
        if(!strcmp(ncandfile.array[i], rxnid))
                return true;
    return false;
}

#define READPARAM(id, type) \
    if (fscanf(fp, #id"=%"#type"\n", &mcprob->id) !=1){ \
        fprintf (stderr, "error: parsing parameter'%s'\n", #id); \
	exit(-1); \
    } // print parameter here if interested

void
loadparams(MCprob *mcprob, const char *filepath){
    FILE *fp = NULL;

    /* open file for reading   */
    if (!(fp = fopen (filepath, "r"))) {
        fprintf (stderr, "error: file open failed '%s'\n", filepath);
	exit(-1);
    }
    READPARAM(objective,s);
    READPARAM(alpha,d);
    READPARAM(beta,d);
    READPARAM(popsize,d);
    READPARAM(ngen,d);
    READPARAM(ncores,d);
    READPARAM(seed,d);
    if (fp) fclose (fp);
}

population
loadpopulation(MCprob *mcprob, const char *inpoppath){
    // If inpoppath is NULL generate a pop randomly

    // Deal with discrepancies between specified popsize and file popsize. Or discrepancies between design variables (e.g. b=0 vs b=1?)

    // Calculate pop objective values?
}

void
writepopulation(MCprob *mcprob, population *outpop, const char *outpoppath){
}

// population runmoea(MCprob *mcprob) // define on separate file

int
main (int argc, char **argv) {
    //TODO: Print only critical info to stdout, print detailed info (including whatever is printed to stdout) to an automatically generated log file. (Or maybe just print everything to log file?)

/* Parse  input */
if (argc < 3){
	fprintf (stderr, "error: insufficient input, usage: %s <problem directory> <parameters> <initial population (optional)>\n", argv[0]);
        return(1);
    }

MCprob mcprob = readproblem(argv[1]);
loadparams(&mcprob, argv[2]);

/* Intialize population */
if (argc == 4){
    char pop_path[256];
    printf("Initial population path: %s \n", argv[3]);
    strcpy(pop_path, argv[3]);
    //inipop = loadpopulation(&mcprob, argv[1]);
}
else{
    printf("Initial population not specified\n");
    //inipop = loadpopulation(&mcprob, NULL);
}

/* Run */
// population pop = runmoea(&mctype, &inipop);
// writepop(mctype, pop);

return(0);
}

