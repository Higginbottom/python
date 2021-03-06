/* lineio.c */
int get_line(FILE *fptr, char line[]);
/* rdpar.c */
int opar(char filename[]);
int add_par(char filename[]);
int cpar(char filename[]);
int rdpar_init(void);
int string_process(char question[], char dummy[]);
int string_process_from_command_line(char question[], char dummy[]);
int string_process_from_file(char question[], char dummy[]);
int rdpar_store_record(char *name, char *value);
int rdpar_save(FILE *file_ptr);
int rdpar_comment(char *format, ...);
int message(char string[]);
int rdstr(char question[], char answer[]);
int rdchar(char question[], char *answer);
int rdint(char question[], int *answer);
int rdint2(char question[], int *answer1, int *answer2);
int rdflo(char question[], float *answer);
int rddoub(char question[], double *answer);
int rdline(char question[], char answer[]);
int get_root(char root[], char total[]);
int rdpar_set_mpi_rank(int rank);
int rdpar_set_verbose(int vlevel);
int rd_extra(char firstword[], double *answer, int *wordlength);
/* log.c */
int Log_init(char *filename);
int Log_append(char *filename);
int Log_close(void);
int Log_set_verbosity(int vlevel);
int Log_print_max(int print_max);
int Log_quit_after_n_errors(int n);
int Log(char *format, ...);
int Log_silent(char *format, ...);
int Error(char *format, ...);
int Error_silent(char *format, ...);
int Shout(char *format, ...);
int sane_check(double x);
int error_count(char *format);
int error_summary(char *message);
int Log_flush(void);
int Log_set_mpi_rank(int rank, int n_mpi);
int Log_parallel(char *format, ...);
int Debug(char *format, ...);
/* synonyms.c */
int check_synonyms(char new_question[], char old_question[]);
