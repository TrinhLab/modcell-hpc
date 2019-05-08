/* Main file, essentially handles IO */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "mcutils.h"
#include "modcell.h"

MCproblem read_problem(const char *problem_dir);
void load_parameters(MCproblem *mcp, const char *filepath);
bool is_not_candidate(Charlist *ncandfile, const char *rxnid);
void write_population(MCproblem *mcp, Population *pop, const char *out_population_path);
void read_population(MCproblem *mcp, Population *pop, const char *population_path);
int get_rxn_idx(MCproblem *mcp, const char *rxn_id);
int get_model_idx(MCproblem *mcp, const char *model_id);

MCproblem
read_problem(const char *problem_dir_path){
    int i,k,j;
    unsigned int n_models = 0;
    unsigned int n_vars;
    char cand_path[256];
    Charlist cand_file, pd;
    MCproblem mcp;
    LPproblem *lp;

    pd = read_dir(problem_dir_path);
    /* Determine number of models and candidates*/
    for (i=0; i < pd.n; i++){
        if(is_extension(pd.array[i], "mps"))
    	    n_models++;
        if(strcmp(pd.array[i], "cand") == 0){
            set_full_path(cand_path, problem_dir_path, pd.array[i]);
    	    cand_file = read_file(cand_path);
        }
    }

    allocate_MCproblem(&mcp, n_models, cand_file.n);

    /* Read deletion candidate names */
    for (j=0; j < mcp.n_vars; j++)
        mcp.individual2id[j] = strdup(cand_file.array[j]);
    free_charlist(cand_file);

    /* Load problems */
    glp_smcp param;
    glp_init_smcp(&param);
    param.msg_lev = GLP_MSG_ERR;
    k=0;
    for (i=0; i < pd.n; i++){
        if(is_extension(pd.array[i], "mps")){
            char *model_name = remove_file_extension(pd.array[i]);
	    mcp.model_names[k] = strdup(model_name);

            char prob_path[256];
            set_full_path(prob_path, problem_dir_path, pd.array[i]);

            lp = &(mcp.lps[k]);
            glp_read_mps(lp->P, GLP_MPS_FILE, NULL, prob_path); // TODO: Silence this or redirect to log:
            /* Calculate basis and solve so subsequent computations are warm-started */
            glp_adv_basis(lp->P, 0);
            glp_simplex(lp->P, &param);
            k++;
        }
    }
    free_charlist(pd);

    /* Gather info to modify LP problems by individuals */
    int col_idx;
    for (k=0; k < mcp.n_models; k++){
        lp = &(mcp.lps[k]);
        char model_path[256];
    	strcpy(model_path, problem_dir_path);
    	strcat(strcat(model_path, glp_get_prob_name(lp->P)), ".ncand");
    	Charlist ncandfile = read_file(model_path);

      	glp_create_index(lp->P);
        for(j=0; j < mcp.n_vars; j++){
            col_idx = glp_find_col(lp->P, mcp.individual2id[j]); /* If glp fails to find the column it will error */
            /* Create maps between indvidual and bounds */
            if(is_not_candidate(&ncandfile, mcp.individual2id[j]))
                lp->cand_col_idx[j] = NOT_CANDIDATE;
            else
                lp->cand_col_idx[j] = col_idx;
            /* Bound info*/
	    lp->cand_col_type[j] = glp_get_col_type(lp->P, col_idx);
	    lp->cand_og_lb[j] = glp_get_col_lb(lp->P, col_idx);
	    lp->cand_og_ub[j] = glp_get_col_ub(lp->P, col_idx);
        }
      	glp_delete_index(lp->P);
    	free_charlist(ncandfile);
    }

return mcp;
}

bool
is_not_candidate(Charlist *ncandfile, const char *rxnid){
    for (int i=0; i < ncandfile->n; i++)
        if(!strcmp(ncandfile->array[i], rxnid))
                return true;
    return false;
}

#define READPARAM(id, type) \
    if (fscanf(fp, #id"=%"#type"\n", &mcp->id) !=1){ \
        fprintf (stderr, "error: parsing parameter'%s'\n", #id); \
        if (fp) fclose (fp); \
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


void
write_population(MCproblem *mcp, Population *pop, const char *out_population_path){
int i,j,k;
FILE *f;
Individual *indv;

if (!(f = fopen (out_population_path, "w"))) {
    fprintf (stderr, "error: file open failed '%s'.", out_population_path);
    exit(-1);
}

fprintf(f, "#METADATA\n");
fprintf(f, "population_size=%d\n", mcp->population_size);
fprintf(f, "alpha=%d\n", mcp->alpha);
fprintf(f, "beta=%d\n", mcp->beta);

for (i=0; i < pop->size; i++){

    indv = &(pop->indv[i]);
    fprintf(f, "#INDIVIDUAL\n");

    fprintf(f, "#DELETIONS\n");
    for (j=0; j < mcp->n_vars; j++)
        if ( indv->deletions[j] == 0 )
            fprintf(f, "%s\n", mcp->individual2id[j]);

    fprintf(f, "#MODULES\n");
    for (k=0; k < mcp->n_models; k++){
       fprintf(f, "%s", mcp->model_names[k]);
        for (j=0; j < mcp->n_vars; j++)
            if (indv->modules[k][j] == 1)
                fprintf(f, ",%s", mcp->individual2id[j]);
        fprintf(f, "\n");
    }

    fprintf(f, "#OBJECTIVES\n");
    for (k=0; k < mcp->n_models; k++)
       fprintf(f, "%s,%f\n", mcp->model_names[k], indv->objectives[k]);
}

fprintf(f, "#ENDFILE\n");
fclose(f);
printf("Output population written to: %s\n",out_population_path);
}

/* Notes: This can be replaced by a hash table, but the look up is only done during IO */
int
get_rxn_idx(MCproblem *mcp, const char *rxn_id){
    for (int j=0; j < mcp->n_vars; j++)
        if(strcmp(mcp->individual2id[j], rxn_id) == 0)
            return j;
    fprintf (stderr, "error: Reaction ID not found: '%s'\n", rxn_id);
    exit(-1);
}

int
get_model_idx(MCproblem *mcp, const char *model_id){
    for (int k=0; k < mcp->n_models; k++)
        if(strcmp(mcp->model_names[k], model_id) == 0)
            return k;
    fprintf (stderr, "error: Model ID not found: '%s'\n", model_id);
    exit(-1);
}

#define READMETADAT(id, type) \
    fscanf(fp, "%s\n", buff); lc++; \
    if (sscanf(buff, #id"=%"#type, &id) !=1){ \
        fprintf (stderr, "error: reading population file parameter: '%s'\n", #id); \
        if (fp) fclose (fp); \
	exit(-1); \
    }

void
read_population(MCproblem *mcp, Population *pop, const char *population_path){
    int population_size, alpha, beta;
    char buff[1000]; //TODO: Is there a way to detect overflow of this buffer?
    int lc = 0;
    bool in_deletions = 0;
    bool in_modules = 0;
    int indv_idx = -1;
    int rxn_idx, model_idx;
    Individual *indv;
    char *token, *string, *tofree;

    FILE *fp = NULL;
    if (!(fp = fopen (population_path, "r"))) {
            fprintf (stderr, "error: file open failed '%s'.", population_path);
            exit(-1);
    }

    /* Read metadata */
    fscanf(fp, "%s\n", buff); lc++;
    assert(strcmp(buff, "#METADATA") == 0);
    READMETADAT(population_size, d)
    READMETADAT(alpha, d)
    READMETADAT(beta, d)

    /* Deal with discrepancies between imput parameters and population file */
    if (population_size > mcp->population_size)
        printf("Population size in parameters (%i) below that of input population (%i), last individuals will be discarded\n", population_size, mcp->population_size);

    if (population_size < mcp->population_size)
        printf("Population size in parameters (%i) above that of input population (%i), missing individuals will be created randomly\n", population_size, mcp->population_size);

    /* Currently ignoring discrepancies in alpha and beta parameters */

    /* Read individuals  */
    while (fscanf(fp, "%s\n", buff) != EOF){
        //printf("%i: %s\n",lc,  buff);
        lc++;

        if(strcmp(buff, "#INDIVIDUAL") == 0){
            indv_idx++;
            if(indv_idx == mcp->population_size)
                break;
            indv = &(pop->indv[indv_idx]);
            set_blank_individual(mcp, indv);
            continue;
        }
        if(strcmp(buff, "#DELETIONS") == 0){
            in_deletions = 1;
            continue;
        }
        if(strcmp(buff, "#MODULES") == 0){
            in_deletions = 0;
            in_modules = 1;
            continue;
        }
        if(strcmp(buff, "#OBJECTIVES") == 0){
            in_modules = 0;
            continue;
        }

        if (in_deletions){
            rxn_idx = get_rxn_idx(mcp, buff);
            indv->deletions[rxn_idx] = 0;
        }

        if (in_modules & (mcp->beta > 0)){
            tofree = string = strdup(buff);
            assert(string != NULL);
            token = strsep(&string, ",");
            model_idx = get_model_idx(mcp, token);
            while ((token = strsep(&string, ",")) != NULL){
                rxn_idx = get_rxn_idx(mcp, token);
                indv->modules[model_idx][rxn_idx] = 1;
            }
        }
    }
    if (fp) fclose(fp);
    if(tofree) free(tofree);

    /* Add extra individuals if needed */
    indv_idx++;
    while(indv_idx < mcp->population_size){
        set_random_individual(mcp, &(pop->indv[indv_idx]));
        indv_idx++;
    }

    /* Calculate objectives */
    for (int i =0; i < mcp->population_size; i++)
            calculate_objectives(mcp, &(pop->indv[i]));
}


int
main (int argc, char **argv){
/* Parse  input */
if ( (argc < 4) | (argc > 5) ){
	fprintf (stderr, "error: insufficient input"
                ", usage: %s <problem_directory> <parameters> <output_population> <initial_population (optional)>\n"
    "Program command-line arguments correspond primarily to file and directory paths:\n"
    "\t- problem_directory: Path to problem directory with a .mps file for each production network, cand file (with a list of reaction candidates) and .ncand file for each production network with a list of fixed reactions per network.\n"
    "\t- parameters: Path to parameter file.\n"
    "\t- output_population: Path to output population file in .pop format. The file does not need to exist (In general it should not exist, since it would be overwritten) but the directory containing the file must exist!\n"
    "\t- input_population (optional): Path to input population file in .pop format. If not included the first population will be initialized randomly.\n"
    "NOTE: Directory for output population must exist.\n"
                ,argv[0]);
        return(1);
    }

MCproblem mcp = read_problem(argv[1]);
load_parameters(&mcp, argv[2]);
printf("--------------------------------------------------\n");

/* Seed global RNG */
pcg32_srandom(mcp.seed, 54u);

/* Intialize population */
Population *initial_population = malloc(sizeof(Population));
allocate_population(&mcp, initial_population, mcp.population_size);

if (argc == 5){
    char pop_path[256];
    printf("Initial population path: %s\n", argv[4]);
    strcpy(pop_path, argv[4]);
    read_population(&mcp, initial_population, pop_path);
}
else{
    printf("Initial population not specified\n");
    set_random_population(&mcp, initial_population);
}

/* Run */
Population output_population = run_moea(&mcp, initial_population);

write_population(&mcp, &output_population, argv[3]);

free_population(&mcp, initial_population);
//free_population(&mcp, &output_population); //FIXME: Uncomment when run_moea is written
return(0);
}

