/* Main file, essentially handles IO */
#pragma GCC diagnostic ignored "-Wmissing-field-initializers" /* Remove warnings about argp structures */
#include <string.h>
#include <assert.h>
#include "mcutils.h"
#include "modcell.h"
#include <stdlib.h>
#include <argp.h>
#include <error.h>
#if MIN_LOG
    #include <unistd.h>
    #include <sys/types.h>
    #include <sys/stat.h>
    #include <fcntl.h>
#endif


/* Macro and function declarations */
#define OPENFILER(filepath)\
    FILE *fp = NULL;\
    if (!(fp = fopen (filepath, "r"))) {\
        fprintf (stderr, "error: file open failed '%s'\n", filepath);\
        exit(-1);\
    }

#define CLOSEFILE if (fp) fclose (fp);

#define READPARAM(id, type, sname) \
    if (fscanf(fp, #id"=%"#type"\n", &(sname->id)) !=1){ \
        fprintf (stderr, "error: parsing parameter'%s'\n", #id); \
        if (fp) fclose (fp); \
        exit(-1); \
    }


#define READMETADAT(id, type) \
    fscanf(fp, "%s\n", buff); lc++; \
    if (sscanf(buff, #id"=%"#type, &id) !=1){ \
        fprintf (stderr, "error: reading population file parameter: '%s'\n", #id); \
        if (fp) fclose (fp); \
        exit(-1); \
    }

MCproblem read_problem(const char *problem_dir);
bool is_not_candidate(Charlist *ncandfile, const char *rxnid);
void write_population(MCproblem *mcp, Population *pop, char *out_population_path);
void read_population(MCproblem *mcp, Population *pop, const char *population_path);
int get_rxn_idx(MCproblem *mcp, const char *rxn_id);
int get_model_idx(MCproblem *mcp, const char *model_id);

extern glp_smcp param;
extern int mpi_pe, mpi_comm_size;

/* Function definitions */

/* CLI */
const char *argp_program_version = MODCELL_V_STRING; /* This variable is set by compiler */
const char *argp_program_bug_address = "https://github.com/trinhlab/modcell-hpc/issues";

/* Program documentation. */
static char doc[] ="Neccessary arguments:\n\t- PROBLEM_DIR: Path to problem directory with a .mps file for each production network, cand file (with a list of reaction candidates) and .ncand file for each production network with a list of fixed reactions per network.\n\t- OUTPUT_FILE: Path to output population file in .pop format. The file does not need to exist (In general it should not exist, since it would be overwritten) but the DIRECTORY containing the file MUST EXIST!\n\
\v modcell-hpc does not check argument bounds, so make sure the values you enter are valid.";

/* A description of the arguments we accept. */
static char args_doc[] = "PROBLEM_DIR OUTPUT_FILE";

/* Keys for options without short-options. */
#define OPT_MINIMIZE_MR  1            /* --minimize_mr */

/* The options we understand. */
static struct argp_option options[] = {
  {"quiet",                     'q', 0,       0, "Don't produce any output" },
  {"initial_population",        'i', "FILE",  0, "Path to input population file in .pop format. If not included the first population will be initialized randomly" },
  {"objective_type",            'd', "STRING",    0, "Design objective. Current options are \"wgcp\"" },
  {"alpha",                     'a', "INT",       0, "Max. number of deletions" },
  {"beta",                      'b', "INT",       0, "Max. number of module reactions" },
  {"seed",                      'r', "INT",       0, "RNG seed. Note that actual seed will be seed + PE_number" },
  {"population_size",           's', "INT",       0, "Number of individuals" },
  {"crossover_probability",     'c', "FLOAT",       0, "Value between 0 and 1 that indicates the chances of crossover for each individual" },
  {"mutation_probability",      'm', "FLOAT",       0, "Value between 0 and 1 that indicates the chances of mutation for each individual" },
  {"migration_interval",        'g', "INT",       0, "Number of generations in between migrations. Note that since migration is asynchronous, this value indicates when migration will be initiated (assuming no migrations are pending, in such case new migration attempts will not be made until ongoing migrations are completed)" },
  {"migration_fraction",        'z', "FLOAT",     0, "Value between 0 and 1. Fraction of the population that will be transfered during migration" },
  {"migration_topology",        'y', "INT",       0, "0: Ring topology, islands communicate as a directed ring graph; 1: Random topology, each migration will send and receive from a random island other than itself." },
  {"migration_policy",          'p', "INT",       0, "0: replace_bottom, the top individuals are sent and the bottom replaced, 1, :replace_sent, the top individuals are sent and replaced; 2, random, Random individuals are sent and replaced. Option 0 maintains the sent individuals in the original population, 1 or 2 do not." },
  {"max_run_time",              't', "INT",       0, "Wall-clock run time in seconds for the main MOEA loop (allow some extra time for IO)" },
  {"n_generations",             'n', "INT",       0, "Maximum number of generations" },
  {"minimize_modules",               OPT_MINIMIZE_MR ,0, 0, "Run module reaction minimizer instead of MOEA"},
  { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char *args[2];     /* arg1 and arg2 */
  char *objective_type, *initial_population;
  int alpha, beta, seed, max_run_time, migration_interval, population_size, verbose, n_generations, migration_policy, migration_topology, minimize_modules;
  float crossover_probability, mutation_probability, migration_fraction;
};

void load_parameters(MCproblem *mcp, struct arguments *arguments);

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'q':
      arguments->verbose = 0;
      break;
    case 'i':
      arguments->initial_population = arg;
      break;
    case 'd':
      arguments->objective_type = arg;
      break;
    case 'a':
      arguments->alpha = atoi(arg);
      break;
    case 'b':
      arguments->beta = atoi(arg);
      break;
    case 'r':
      arguments->seed = atoi(arg);
      break;
    case 's':
      arguments->population_size = atoi(arg);
      break;
    case 'c':
      arguments->crossover_probability = atof(arg);
      break;
    case 'm':
      arguments->mutation_probability = atof(arg);
      break;
    case 'g':
      arguments->migration_interval = atoi(arg);
      break;
    case 'z':
      arguments->migration_fraction = atof(arg);
      break;
    case 'y':
      arguments->migration_topology = atoi(arg);
      break;
    case 'p':
      arguments->migration_policy = atoi(arg);
      break;
    case 't':
      arguments->max_run_time = atoi(arg);
      break;
    case 'n':
      arguments->n_generations = atoi(arg);
      break;
    case OPT_MINIMIZE_MR:
      arguments->minimize_modules = 1;
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num >= 2) /* Too many arguments. */
        argp_usage (state);
      arguments->args[state->arg_num] = arg;
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2) /* Not enough arguments. */
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

void
load_parameters(MCproblem *mcp, struct arguments *arguments)
{
    mcp->verbose = arguments->verbose;
    strcpy(mcp->objective_type, arguments->objective_type);
    mcp->alpha = arguments->alpha;
    mcp->beta = arguments->beta;
    mcp->seed = arguments->seed;
    mcp->population_size = arguments->population_size;
    mcp->crossover_probability = arguments->crossover_probability;
    mcp->mutation_probability = arguments->mutation_probability;
    mcp->migration_interval = arguments->migration_interval;
    mcp->migration_size = (int)(arguments->migration_fraction * arguments->population_size);
    mcp->max_run_time = arguments->max_run_time;
    mcp->n_generations = arguments->n_generations;
    mcp->migration_topology = arguments->migration_topology;
    mcp->migration_policy = arguments->migration_policy;
    /* Indicate if module reactions are used */
    mcp->use_modules = arguments->beta > 0;
}

/* CLI done */

MCproblem
read_problem(const char *problem_dir_path)
{
    int i,k,j,col_idx;
    unsigned int n_models = 0;
    char cand_path[256], prob_path[256], model_path[256], mparam_path[256], prod_rxn_name[256];
    Charlist cand_file={NULL}, pd;
    MCproblem mcp;
    LPproblem *lp;

    pd = read_dir(problem_dir_path);

    /* Determine number of models and candidates*/
    for (i=0; i < pd.n; i++) {
        if(is_extension(pd.array[i], "mps"))
            n_models++;
        if(strcmp(pd.array[i], "cand") == 0) {
            set_full_path(cand_path, problem_dir_path, pd.array[i]);
            cand_file = read_file(cand_path);
        }
    }

    allocate_MCproblem(&mcp, n_models, cand_file.n);

    /* Read deletion candidate names */
    for (j=0; j < mcp.n_vars; j++)
        mcp.individual2id[j] = strdup(cand_file.array[j]);
    free_charlist(cand_file);

    #if MIN_LOG
        /* Temporary redirection of stdout to silence glpk */
        int bak, new;
        fflush(stdout);
        bak = dup(1);
        new = open("/dev/null", O_WRONLY);
        dup2(new, 1);
        close(new);
    #endif

    /* Load problems */
    k=0;
    for (i=0; i < pd.n; i++) {
        if(is_extension(pd.array[i], "mps")) {
            char *model_name = remove_file_extension(pd.array[i]);
            mcp.model_names[k] = strdup(model_name);

            set_full_path(prob_path, problem_dir_path, pd.array[i]);

            lp = &(mcp.lps[k]);
            glp_read_mps(lp->P, GLP_MPS_FILE, NULL, prob_path);
            /* Calculate basis and solve so subsequent computations are warm-started */
            glp_adv_basis(lp->P, 0);
            glp_simplex(lp->P, &param);
            k++;
        }
    }
    free_charlist(pd);

    #if MIN_LOG
        fflush(stdout);
        dup2(bak, 1);
        close(bak);
    #endif

    /* Gather info to modify LP problems by individuals */
    for (k=0; k < mcp.n_models; k++){
        lp = &(mcp.lps[k]);
        strcpy(model_path, problem_dir_path);
        strcat(strcat(model_path, glp_get_prob_name(lp->P)), ".ncand");
        Charlist ncandfile = read_file(model_path);

        glp_create_index(lp->P);
        for(j=0; j < mcp.n_vars; j++) {
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

        /* Other parameters */
        strcpy(mparam_path, problem_dir_path);
        strcat(strcat(mparam_path, glp_get_prob_name(lp->P)), ".param");

        OPENFILER(mparam_path)
            if (fscanf(fp, "prod_rxn_name=%s\n", prod_rxn_name) !=1) {
                fprintf (stderr, "error: parsing parameter prod_rxn_name\n"); exit(-1);}

        lp->prod_col_idx = glp_find_col(lp->P, prod_rxn_name);
        READPARAM(max_prod_growth,lf,lp)
            CLOSEFILE

            /* Cleanup */
            glp_delete_index(lp->P);
        free_charlist(ncandfile);

        /* Objective values without deletions */
        glp_simplex(lp->P, &param);
        lp->no_deletion_objective = glp_get_col_prim(lp->P, lp->prod_col_idx)/lp->max_prod_growth; //FIXME: Assumes wGCP. Needs to be calculated after parameters are parsed.
    }

    return mcp;
}

bool
is_not_candidate(Charlist *ncandfile, const char *rxnid)
{
    for (int i=0; i < ncandfile->n; i++)
        if(!strcmp(ncandfile->array[i], rxnid))
            return true;
    return false;
}



void
write_population(MCproblem *mcp, Population *pop, char *out_population_path)
{
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

    for (i=0; i < pop->size; i++) {
        indv = &(pop->indv[i]);
        fprintf(f, "#INDIVIDUAL\n");

        fprintf(f, "#DELETIONS\n");
        for (j=0; j < mcp->n_vars; j++)
            if ( indv->deletions[j] == DELETED_RXN )
                fprintf(f, "%s\n", mcp->individual2id[j]);

        fprintf(f, "#MODULES\n");
        for (k=0; k < mcp->n_models; k++) {
            fprintf(f, "%s", mcp->model_names[k]);
            for (j=0; j < mcp->n_vars; j++)
                if (mcp->use_modules)
                    if (indv->modules[k*mcp->n_vars + j] == MODULE_RXN)
                        fprintf(f, ",%s", mcp->individual2id[j]);
            fprintf(f, "\n");
        }

        fprintf(f, "#OBJECTIVES\n");
        for (k=0; k < mcp->n_models; k++)
            fprintf(f, "%s,%lf\n", mcp->model_names[k], indv->objectives[k]);
    }

    fprintf(f, "#ENDFILE\n");
    fclose(f);
    printf("Output population written to: %s\n",out_population_path);
}

/* Notes: This can be replaced by a hash table, but the look up is only done during IO */
int
get_rxn_idx(MCproblem *mcp, const char *rxn_id)
{
    for (int j=0; j < mcp->n_vars; j++)
        if(strcmp(mcp->individual2id[j], rxn_id) == 0)
            return j;
    fprintf (stderr, "error: Reaction ID not found: '%s'\n", rxn_id);
    exit(-1);
}

int
get_model_idx(MCproblem *mcp, const char *model_id)
{
    for (int k=0; k < mcp->n_models; k++)
        if(strcmp(mcp->model_names[k], model_id) == 0)
            return k;
    fprintf (stderr, "error: Model ID not found: '%s'\n", model_id);
    exit(-1);
}


/* Loads a population file:
 * Notes:
 *      - Objectives could be calculated here and checked for consitency, currently objectives are not calculated until the moea procdure.
 */
void
read_population(MCproblem *mcp, Population *pop, const char *population_path)
{
    int population_size, alpha, beta;
    char buff[1000]; //TODO: Is there a way to detect overflow of this buffer?
    int lc = 0;
    bool in_deletions = 0, in_modules = 0;
    int indv_idx = -1, rxn_idx, model_idx;
    Individual *indv ={NULL};
    char *token, *string, *tofree=NULL;

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
        printf("Population size in parameters (%i) below that of input population (%i), last individuals will be discarded\n", mcp->population_size, population_size);

    if (population_size < mcp->population_size)
        printf("Population size in parameters (%i) above that of input population (%i), missing individuals will be created randomly\n", mcp->population_size,  population_size);

    /* Currently ignoring discrepancies in alpha and beta parameters */

    /* Read individuals  */
    while (fscanf(fp, "%s\n", buff) != EOF) {
        lc++;

        if(strcmp(buff, "#INDIVIDUAL") == 0) {
            indv_idx++;
            if(indv_idx == mcp->population_size)
                break;
            indv = &(pop->indv[indv_idx]);
            set_blank_individual(mcp, indv);
            continue;
        }
        if(strcmp(buff, "#DELETIONS") == 0) {
            in_deletions = 1;
            continue;
        }
        if(strcmp(buff, "#MODULES") == 0) {
            in_deletions = 0;
            in_modules = 1;
            continue;
        }
        if(strcmp(buff, "#OBJECTIVES") == 0) {
            in_modules = 0;
            continue;
        }

        if (in_deletions) {
            rxn_idx = get_rxn_idx(mcp, buff);
            indv->deletions[rxn_idx] = DELETED_RXN;
        }

        if (in_modules && (mcp->use_modules)) {
            tofree = string = strdup(buff);
            assert(string != NULL);
            token = strsep(&string, ",");
            model_idx = get_model_idx(mcp, token);
            while ((token = strsep(&string, ",")) != NULL){
                rxn_idx = get_rxn_idx(mcp, token);
                indv->modules[model_idx*mcp->n_vars + rxn_idx] = MODULE_RXN;
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

}


int
main (int argc, char **argv)
{
    printf("Modcell version: %s\n", MODCELL_V_STRING);

    /* Parse  input */
    struct arguments arguments;

    /* Default values. */
    arguments.verbose = 1;
    arguments.initial_population = "";
    arguments.objective_type = "wgcp";
    arguments.alpha = 5;
    arguments.beta = 0;
    arguments.seed = 0;
    arguments.population_size = 200;
    arguments.crossover_probability = 0.8;
    arguments.mutation_probability = 0.05;
    arguments.migration_interval = 20;
    arguments.migration_fraction = 0.1;
    arguments.max_run_time = 10000000;
    arguments.n_generations = 500;
    arguments.migration_policy = 0;
    arguments.migration_topology = 0;
    arguments.minimize_modules = 0;

    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    /* Intialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_pe);

    if (mpi_pe == 0) printf("Comm size: %d.\n", mpi_comm_size);

    /* Intialize global GLPK parameters*/
    glp_init_smcp(&param);
    param.msg_lev = LP_MSG_LEV;
    param.tm_lim = LP_TIME_LIMIT_MILISEC;

    MCproblem mcp = read_problem(arguments.args[0]);
    load_parameters(&mcp, &arguments);
    fflush(stdout);

    /* Seed global RNG */
    pcg32_srandom(mcp.seed+mpi_pe, 54u);

    /* Intialize population */
    Population *initial_population = malloc(sizeof(Population));
    allocate_population(&mcp, initial_population, mcp.population_size);
    if (arguments.initial_population[0] == '\0') {
        if (mpi_pe == 0)  printf("(PE=0) Initial population not specified (initialize at random)\n");
        set_random_population(&mcp, initial_population);
    }
    else {
        if (mpi_pe == 0) printf("(PE=0) Reading initial population (path: %s)...", arguments.initial_population);
        read_population(&mcp, initial_population, arguments.initial_population);
        printf("done\n");
    }

    /* Run */
    if (arguments.minimize_modules)  {
        printf("Performing module minimization. MOEA will NOT run.\n");
            if (mpi_comm_size > 1) {
                fprintf (stderr, "error: Module minimization does not run in parallel so invoke without MPI.\n");
                exit(-1);
            }
        minimize_mr(&mcp, initial_population);
    }
    else
        run_moea(&mcp, initial_population);

    /* Write ouput */
    char pop_path[256];
    if (mpi_comm_size > 1)
        sprintf(pop_path, "%s_%i", arguments.args[1], mpi_pe);
    else
        sprintf(pop_path, "%s", arguments.args[1]);
    write_population(&mcp, initial_population, pop_path);

    /* Wait for all processes before exiting (Avoid attempts to communicate with finished processes).*/
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_pe == 0) printf("Barrier reached, all processes exiting...\n");

    /* Cleanup */
    free_population(&mcp, initial_population);
    MPI_Finalize();

    return(0);
}

