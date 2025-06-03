#ifndef CMDLINE_H
#define CMDLINE_H
#include <string>

#define CParty_CMDLINE_PACKAGE_NAME "CParty"

#define CParty_CMDLINE_VERSION "1.0"


// The restricted structure
extern std::string input_struct;

// The input file
extern std::string input_file;

// The output file
extern std::string output_file;

// Number of Suboptimal structures to print
extern int subopt;

// Number for dangle model
extern int dangle_model;

// The parameter file location
extern std::string parameter_file;

// Number of Suboptimal structures to print
extern int samples;
// The shape file
// extern std::string shape_file;



/** @brief Where the command line options are stored */
struct args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  const char *input_structure_help; /**< @brief Give restricted structure as input help description.  */
  const char *input_file_help; /**< @brief Give input file as input help description.  */
  const char *output_file_help; /**< @brief Give output file as input help description.  */
  const char *subopt_help; /**< @brief Give a number of suboptimals to print  */
  const char *pk_free_help; /**< @brief Specify if the structures calculated should be pseudoknot-free  */
  const char *pk_only_help; /**< @brief Specify if the structures added should all cross the input   */
  const char *dangles_help; /**< @brief Specify the dangle model*/
  const char *paramFile_help; /**< @brief Use a separate parameter list */
  const char *samples_help; /**< @brief Specify the number of samples for the stochastic backtracking (default 1000).  */
  // const char *shape_help; /**< @brief Give shape file as additional input help description.  */
  const char *noConv_help; /**< @brief Turn off automated conversion to RNA help description.  */


  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_structure_given ;	/**< @brief Whether restricted structure was given.  */
  unsigned int input_file_given ;	/**< @brief Whether restricted structure was given.  */
  unsigned int output_file_given ;	/**< @brief Whether restricted structure was given.  */
  unsigned int subopt_given ;	/**< @brief Whether suboptimals was given.  */
  unsigned int pk_free_given ;	/**< @brief Whether pk_free was given.  */
  unsigned int pk_only_given ;	/**< @brief Whether pk_only was given.  */
  unsigned int dangles_given ;  /**< @brief Whether dangle model was given.  */
  unsigned int paramFile_given ; /** <@brief whether a parameter file was given */
  unsigned int samples_given ; /**< @brief Whether samples was given.  */
  // unsigned int shape_given ; /**< @brief Whether shape was given.  */
  unsigned int noConv_given ;	/**< @brief Whether noConv was given.  */

  char **inputs ; /**< @brief unnamed options (options without names) */
  unsigned inputs_num ; /**< @brief unnamed options number */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];


/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,struct args_info *args_info);


/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv, struct args_info *args_info);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct args_info *args_info);

/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct args_info *args_info);






#endif