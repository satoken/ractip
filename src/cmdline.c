/*
  File autogenerated by gengetopt version 2.22.6
  generated with the following command:
  gengetopt -u

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "RactIP: RNA-RNA interation prediction using integer programming.";

const char *gengetopt_args_info_usage = "Usage: " CMDLINE_PARSER_PACKAGE " [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_versiontext = "";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_full_help[] = {
  "  -h, --help                 Print help and exit",
  "      --full-help            Print help, including hidden options, and exit",
  "  -V, --version              Print version and exit",
  "  -a, --alpha=FLOAT          weight for hybridization  (default=`0.7')",
  "  -b, --beta=FLOAT           weight for accessibility  (default=`0.0')",
  "  -t, --fold-th=FLOAT        Threshold for base-pairing probabilities\n                               (default=`0.5')",
  "  -u, --hybridize-th=FLOAT   Threshold for hybridazation probabilities\n                               (default=`0.1')",
  "  -s, --acc-th=FLOAT         Threshold for accessible probabilities\n                               (default=`0.003')",
  "      --acc-max              optimize for accessibility instead of internal\n                               secondary structures  (default=off)",
  "      --acc-max-ss           additional prediction of interanal secondary\n                               structures  (default=off)",
  "      --acc-num=INT          the number of accessible regions (0=unlimited)\n                               (default=`1')",
  "      --max-w=INT            Maximum length of accessible regions\n                               (default=`15')",
  "      --min-w=INT            Minimum length of accessible regions\n                               (default=`5')",
  "      --zscore=INT           Calculate z-score via dishuffling (0=no shuffling,\n                               1=1st seq only, 2=2nd seq only, or 12=both)\n                               (default=`0')",
  "      --num-shuffling=INT    The number of shuffling  (default=`1000')",
  "      --seed=INT             Seed for random number generator  (default=`0')",
  "  -c, --use-constraint       Use structure constraints  (default=off)",
  "      --force-constraint     Enforce structure constraints  (default=off)",
  "      --allow-isolated       Allow isolated base-pairs  (default=off)",
  "  -e, --show-energy          calculate the free energy of the predicted joint\n                               structure  (default=off)",
  "  -P, --param-file=FILENAME  Read the energy parameter file for Vienna RNA\n                               package",
  "      --no-pk                do not use the constraints for interenal\n                               pseudoknots  (default=off)",
  "  -r, --rip=FILENAME         Import posterior probabilities from the result of\n                               RIP",
  "      --no-bl                do not use BL parameters  (default=off)",
    0
};

static void
init_help_array(void)
{
  gengetopt_args_info_help[0] = gengetopt_args_info_full_help[0];
  gengetopt_args_info_help[1] = gengetopt_args_info_full_help[1];
  gengetopt_args_info_help[2] = gengetopt_args_info_full_help[2];
  gengetopt_args_info_help[3] = gengetopt_args_info_full_help[3];
  gengetopt_args_info_help[4] = gengetopt_args_info_full_help[4];
  gengetopt_args_info_help[5] = gengetopt_args_info_full_help[5];
  gengetopt_args_info_help[6] = gengetopt_args_info_full_help[6];
  gengetopt_args_info_help[7] = gengetopt_args_info_full_help[7];
  gengetopt_args_info_help[8] = gengetopt_args_info_full_help[8];
  gengetopt_args_info_help[9] = gengetopt_args_info_full_help[9];
  gengetopt_args_info_help[10] = gengetopt_args_info_full_help[10];
  gengetopt_args_info_help[11] = gengetopt_args_info_full_help[11];
  gengetopt_args_info_help[12] = gengetopt_args_info_full_help[12];
  gengetopt_args_info_help[13] = gengetopt_args_info_full_help[13];
  gengetopt_args_info_help[14] = gengetopt_args_info_full_help[14];
  gengetopt_args_info_help[15] = gengetopt_args_info_full_help[15];
  gengetopt_args_info_help[16] = gengetopt_args_info_full_help[18];
  gengetopt_args_info_help[17] = gengetopt_args_info_full_help[19];
  gengetopt_args_info_help[18] = gengetopt_args_info_full_help[20];
  gengetopt_args_info_help[19] = 0; 
  
}

const char *gengetopt_args_info_help[20];

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_FLOAT
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->full_help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->alpha_given = 0 ;
  args_info->beta_given = 0 ;
  args_info->fold_th_given = 0 ;
  args_info->hybridize_th_given = 0 ;
  args_info->acc_th_given = 0 ;
  args_info->acc_max_given = 0 ;
  args_info->acc_max_ss_given = 0 ;
  args_info->acc_num_given = 0 ;
  args_info->max_w_given = 0 ;
  args_info->min_w_given = 0 ;
  args_info->zscore_given = 0 ;
  args_info->num_shuffling_given = 0 ;
  args_info->seed_given = 0 ;
  args_info->use_constraint_given = 0 ;
  args_info->force_constraint_given = 0 ;
  args_info->allow_isolated_given = 0 ;
  args_info->show_energy_given = 0 ;
  args_info->param_file_given = 0 ;
  args_info->no_pk_given = 0 ;
  args_info->rip_given = 0 ;
  args_info->no_bl_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->alpha_arg = 0.7;
  args_info->alpha_orig = NULL;
  args_info->beta_arg = 0.0;
  args_info->beta_orig = NULL;
  args_info->fold_th_arg = 0.5;
  args_info->fold_th_orig = NULL;
  args_info->hybridize_th_arg = 0.1;
  args_info->hybridize_th_orig = NULL;
  args_info->acc_th_arg = 0.003;
  args_info->acc_th_orig = NULL;
  args_info->acc_max_flag = 0;
  args_info->acc_max_ss_flag = 0;
  args_info->acc_num_arg = 1;
  args_info->acc_num_orig = NULL;
  args_info->max_w_arg = 15;
  args_info->max_w_orig = NULL;
  args_info->min_w_arg = 5;
  args_info->min_w_orig = NULL;
  args_info->zscore_arg = 0;
  args_info->zscore_orig = NULL;
  args_info->num_shuffling_arg = 1000;
  args_info->num_shuffling_orig = NULL;
  args_info->seed_arg = 0;
  args_info->seed_orig = NULL;
  args_info->use_constraint_flag = 0;
  args_info->force_constraint_flag = 0;
  args_info->allow_isolated_flag = 0;
  args_info->show_energy_flag = 0;
  args_info->param_file_arg = NULL;
  args_info->param_file_orig = NULL;
  args_info->no_pk_flag = 0;
  args_info->rip_arg = NULL;
  args_info->rip_orig = NULL;
  args_info->no_bl_flag = 0;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{

  init_help_array(); 
  args_info->help_help = gengetopt_args_info_full_help[0] ;
  args_info->full_help_help = gengetopt_args_info_full_help[1] ;
  args_info->version_help = gengetopt_args_info_full_help[2] ;
  args_info->alpha_help = gengetopt_args_info_full_help[3] ;
  args_info->beta_help = gengetopt_args_info_full_help[4] ;
  args_info->fold_th_help = gengetopt_args_info_full_help[5] ;
  args_info->hybridize_th_help = gengetopt_args_info_full_help[6] ;
  args_info->acc_th_help = gengetopt_args_info_full_help[7] ;
  args_info->acc_max_help = gengetopt_args_info_full_help[8] ;
  args_info->acc_max_ss_help = gengetopt_args_info_full_help[9] ;
  args_info->acc_num_help = gengetopt_args_info_full_help[10] ;
  args_info->max_w_help = gengetopt_args_info_full_help[11] ;
  args_info->min_w_help = gengetopt_args_info_full_help[12] ;
  args_info->zscore_help = gengetopt_args_info_full_help[13] ;
  args_info->num_shuffling_help = gengetopt_args_info_full_help[14] ;
  args_info->seed_help = gengetopt_args_info_full_help[15] ;
  args_info->use_constraint_help = gengetopt_args_info_full_help[16] ;
  args_info->force_constraint_help = gengetopt_args_info_full_help[17] ;
  args_info->allow_isolated_help = gengetopt_args_info_full_help[18] ;
  args_info->show_energy_help = gengetopt_args_info_full_help[19] ;
  args_info->param_file_help = gengetopt_args_info_full_help[20] ;
  args_info->no_pk_help = gengetopt_args_info_full_help[21] ;
  args_info->rip_help = gengetopt_args_info_full_help[22] ;
  args_info->no_bl_help = gengetopt_args_info_full_help[23] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(CMDLINE_PARSER_PACKAGE_NAME) ? CMDLINE_PARSER_PACKAGE_NAME : CMDLINE_PARSER_PACKAGE),
     CMDLINE_PARSER_VERSION);

  if (strlen(gengetopt_args_info_versiontext) > 0)
    printf("\n%s\n", gengetopt_args_info_versiontext);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n\n", gengetopt_args_info_description);
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_print_full_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_full_help[i])
    printf("%s\n", gengetopt_args_info_full_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void)
{
  struct cmdline_parser_params *params = 
    (struct cmdline_parser_params *)malloc(sizeof(struct cmdline_parser_params));
  cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->alpha_orig));
  free_string_field (&(args_info->beta_orig));
  free_string_field (&(args_info->fold_th_orig));
  free_string_field (&(args_info->hybridize_th_orig));
  free_string_field (&(args_info->acc_th_orig));
  free_string_field (&(args_info->acc_num_orig));
  free_string_field (&(args_info->max_w_orig));
  free_string_field (&(args_info->min_w_orig));
  free_string_field (&(args_info->zscore_orig));
  free_string_field (&(args_info->num_shuffling_orig));
  free_string_field (&(args_info->seed_orig));
  free_string_field (&(args_info->param_file_arg));
  free_string_field (&(args_info->param_file_orig));
  free_string_field (&(args_info->rip_arg));
  free_string_field (&(args_info->rip_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->full_help_given)
    write_into_file(outfile, "full-help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->alpha_given)
    write_into_file(outfile, "alpha", args_info->alpha_orig, 0);
  if (args_info->beta_given)
    write_into_file(outfile, "beta", args_info->beta_orig, 0);
  if (args_info->fold_th_given)
    write_into_file(outfile, "fold-th", args_info->fold_th_orig, 0);
  if (args_info->hybridize_th_given)
    write_into_file(outfile, "hybridize-th", args_info->hybridize_th_orig, 0);
  if (args_info->acc_th_given)
    write_into_file(outfile, "acc-th", args_info->acc_th_orig, 0);
  if (args_info->acc_max_given)
    write_into_file(outfile, "acc-max", 0, 0 );
  if (args_info->acc_max_ss_given)
    write_into_file(outfile, "acc-max-ss", 0, 0 );
  if (args_info->acc_num_given)
    write_into_file(outfile, "acc-num", args_info->acc_num_orig, 0);
  if (args_info->max_w_given)
    write_into_file(outfile, "max-w", args_info->max_w_orig, 0);
  if (args_info->min_w_given)
    write_into_file(outfile, "min-w", args_info->min_w_orig, 0);
  if (args_info->zscore_given)
    write_into_file(outfile, "zscore", args_info->zscore_orig, 0);
  if (args_info->num_shuffling_given)
    write_into_file(outfile, "num-shuffling", args_info->num_shuffling_orig, 0);
  if (args_info->seed_given)
    write_into_file(outfile, "seed", args_info->seed_orig, 0);
  if (args_info->use_constraint_given)
    write_into_file(outfile, "use-constraint", 0, 0 );
  if (args_info->force_constraint_given)
    write_into_file(outfile, "force-constraint", 0, 0 );
  if (args_info->allow_isolated_given)
    write_into_file(outfile, "allow-isolated", 0, 0 );
  if (args_info->show_energy_given)
    write_into_file(outfile, "show-energy", 0, 0 );
  if (args_info->param_file_given)
    write_into_file(outfile, "param-file", args_info->param_file_orig, 0);
  if (args_info->no_pk_given)
    write_into_file(outfile, "no-pk", 0, 0 );
  if (args_info->rip_given)
    write_into_file(outfile, "rip", args_info->rip_orig, 0);
  if (args_info->no_bl_given)
    write_into_file(outfile, "no-bl", 0, 0 );
  

  i = EXIT_SUCCESS;
  return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char **argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char **argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser2 (int argc, char **argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  FIX_UNUSED (args_info);
  FIX_UNUSED (prog_name);
  return EXIT_SUCCESS;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_FLOAT:
    if (val) *((float *)field) = (float)strtod (val, &stop_char);
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_FLOAT:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
cmdline_parser_internal (
  int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error_occurred = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "full-help",	0, NULL, 0 },
        { "version",	0, NULL, 'V' },
        { "alpha",	1, NULL, 'a' },
        { "beta",	1, NULL, 'b' },
        { "fold-th",	1, NULL, 't' },
        { "hybridize-th",	1, NULL, 'u' },
        { "acc-th",	1, NULL, 's' },
        { "acc-max",	0, NULL, 0 },
        { "acc-max-ss",	0, NULL, 0 },
        { "acc-num",	1, NULL, 0 },
        { "max-w",	1, NULL, 0 },
        { "min-w",	1, NULL, 0 },
        { "zscore",	1, NULL, 0 },
        { "num-shuffling",	1, NULL, 0 },
        { "seed",	1, NULL, 0 },
        { "use-constraint",	0, NULL, 'c' },
        { "force-constraint",	0, NULL, 0 },
        { "allow-isolated",	0, NULL, 0 },
        { "show-energy",	0, NULL, 'e' },
        { "param-file",	1, NULL, 'P' },
        { "no-pk",	0, NULL, 0 },
        { "rip",	1, NULL, 'r' },
        { "no-bl",	0, NULL, 0 },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVa:b:t:u:s:ceP:r:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          cmdline_parser_print_version ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'a':	/* weight for hybridization.  */
        
        
          if (update_arg( (void *)&(args_info->alpha_arg), 
               &(args_info->alpha_orig), &(args_info->alpha_given),
              &(local_args_info.alpha_given), optarg, 0, "0.7", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "alpha", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* weight for accessibility.  */
        
        
          if (update_arg( (void *)&(args_info->beta_arg), 
               &(args_info->beta_orig), &(args_info->beta_given),
              &(local_args_info.beta_given), optarg, 0, "0.0", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "beta", 'b',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Threshold for base-pairing probabilities.  */
        
        
          if (update_arg( (void *)&(args_info->fold_th_arg), 
               &(args_info->fold_th_orig), &(args_info->fold_th_given),
              &(local_args_info.fold_th_given), optarg, 0, "0.5", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "fold-th", 't',
              additional_error))
            goto failure;
        
          break;
        case 'u':	/* Threshold for hybridazation probabilities.  */
        
        
          if (update_arg( (void *)&(args_info->hybridize_th_arg), 
               &(args_info->hybridize_th_orig), &(args_info->hybridize_th_given),
              &(local_args_info.hybridize_th_given), optarg, 0, "0.1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "hybridize-th", 'u',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Threshold for accessible probabilities.  */
        
        
          if (update_arg( (void *)&(args_info->acc_th_arg), 
               &(args_info->acc_th_orig), &(args_info->acc_th_given),
              &(local_args_info.acc_th_given), optarg, 0, "0.003", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "acc-th", 's',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Use structure constraints.  */
        
        
          if (update_arg((void *)&(args_info->use_constraint_flag), 0, &(args_info->use_constraint_given),
              &(local_args_info.use_constraint_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "use-constraint", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* calculate the free energy of the predicted joint structure.  */
        
        
          if (update_arg((void *)&(args_info->show_energy_flag), 0, &(args_info->show_energy_given),
              &(local_args_info.show_energy_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "show-energy", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'P':	/* Read the energy parameter file for Vienna RNA package.  */
        
        
          if (update_arg( (void *)&(args_info->param_file_arg), 
               &(args_info->param_file_orig), &(args_info->param_file_given),
              &(local_args_info.param_file_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "param-file", 'P',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Import posterior probabilities from the result of RIP.  */
        
        
          if (update_arg( (void *)&(args_info->rip_arg), 
               &(args_info->rip_orig), &(args_info->rip_given),
              &(local_args_info.rip_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "rip", 'r',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          if (strcmp (long_options[option_index].name, "full-help") == 0) {
            cmdline_parser_print_full_help ();
            cmdline_parser_free (&local_args_info);
            exit (EXIT_SUCCESS);
          }

          /* optimize for accessibility instead of internal secondary structures.  */
          if (strcmp (long_options[option_index].name, "acc-max") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->acc_max_flag), 0, &(args_info->acc_max_given),
                &(local_args_info.acc_max_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "acc-max", '-',
                additional_error))
              goto failure;
          
          }
          /* additional prediction of interanal secondary structures.  */
          else if (strcmp (long_options[option_index].name, "acc-max-ss") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->acc_max_ss_flag), 0, &(args_info->acc_max_ss_given),
                &(local_args_info.acc_max_ss_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "acc-max-ss", '-',
                additional_error))
              goto failure;
          
          }
          /* the number of accessible regions (0=unlimited).  */
          else if (strcmp (long_options[option_index].name, "acc-num") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->acc_num_arg), 
                 &(args_info->acc_num_orig), &(args_info->acc_num_given),
                &(local_args_info.acc_num_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "acc-num", '-',
                additional_error))
              goto failure;
          
          }
          /* Maximum length of accessible regions.  */
          else if (strcmp (long_options[option_index].name, "max-w") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->max_w_arg), 
                 &(args_info->max_w_orig), &(args_info->max_w_given),
                &(local_args_info.max_w_given), optarg, 0, "15", ARG_INT,
                check_ambiguity, override, 0, 0,
                "max-w", '-',
                additional_error))
              goto failure;
          
          }
          /* Minimum length of accessible regions.  */
          else if (strcmp (long_options[option_index].name, "min-w") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->min_w_arg), 
                 &(args_info->min_w_orig), &(args_info->min_w_given),
                &(local_args_info.min_w_given), optarg, 0, "5", ARG_INT,
                check_ambiguity, override, 0, 0,
                "min-w", '-',
                additional_error))
              goto failure;
          
          }
          /* Calculate z-score via dishuffling (0=no shuffling, 1=1st seq only, 2=2nd seq only, or 12=both).  */
          else if (strcmp (long_options[option_index].name, "zscore") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->zscore_arg), 
                 &(args_info->zscore_orig), &(args_info->zscore_given),
                &(local_args_info.zscore_given), optarg, 0, "0", ARG_INT,
                check_ambiguity, override, 0, 0,
                "zscore", '-',
                additional_error))
              goto failure;
          
          }
          /* The number of shuffling.  */
          else if (strcmp (long_options[option_index].name, "num-shuffling") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->num_shuffling_arg), 
                 &(args_info->num_shuffling_orig), &(args_info->num_shuffling_given),
                &(local_args_info.num_shuffling_given), optarg, 0, "1000", ARG_INT,
                check_ambiguity, override, 0, 0,
                "num-shuffling", '-',
                additional_error))
              goto failure;
          
          }
          /* Seed for random number generator.  */
          else if (strcmp (long_options[option_index].name, "seed") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->seed_arg), 
                 &(args_info->seed_orig), &(args_info->seed_given),
                &(local_args_info.seed_given), optarg, 0, "0", ARG_INT,
                check_ambiguity, override, 0, 0,
                "seed", '-',
                additional_error))
              goto failure;
          
          }
          /* Enforce structure constraints.  */
          else if (strcmp (long_options[option_index].name, "force-constraint") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->force_constraint_flag), 0, &(args_info->force_constraint_given),
                &(local_args_info.force_constraint_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "force-constraint", '-',
                additional_error))
              goto failure;
          
          }
          /* Allow isolated base-pairs.  */
          else if (strcmp (long_options[option_index].name, "allow-isolated") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->allow_isolated_flag), 0, &(args_info->allow_isolated_given),
                &(local_args_info.allow_isolated_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "allow-isolated", '-',
                additional_error))
              goto failure;
          
          }
          /* do not use the constraints for interenal pseudoknots.  */
          else if (strcmp (long_options[option_index].name, "no-pk") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->no_pk_flag), 0, &(args_info->no_pk_given),
                &(local_args_info.no_pk_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "no-pk", '-',
                additional_error))
              goto failure;
          
          }
          /* do not use BL parameters.  */
          else if (strcmp (long_options[option_index].name, "no-bl") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->no_bl_flag), 0, &(args_info->no_bl_given),
                &(local_args_info.no_bl_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "no-bl", '-',
                additional_error))
              goto failure;
          
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */




  cmdline_parser_release (&local_args_info);

  if ( error_occurred )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
