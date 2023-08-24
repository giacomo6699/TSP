#include "tsp_solver.h"
#include <time.h>
#include <cplex.h> 

double random01();

void print_uncovered(int *uncovered, int uncovered_len);

void print_succ(int *succ, int nnodes);

void print_travel(int *succ, instance *inst);

void print_instance(instance *inst);

void print_error(const char *err);

double compute_sol_cost(instance *inst, int *succ, int succ_size);

double eucl_dist(int i, int j, instance *inst);

int is_fractional(double x);

int is_all_integer(int n, const double *x);

void double_vector_copy(int n, const double *from, double *to);

void choose_starter(int *starter, instance *inst);

void choose_random_starter(int *starter, instance *inst);

void choose_starters(int *a, int *b, instance *inst);

void choose_random_starters(int *a, int *b, instance *inst);

void reverse_edges(int b, int a_prime, int *succ);

double two_opt_move(instance *inst, int *succ, double curr_cost);

double two_opt_improvement(instance *inst, int *succ, double curr_cost, double timelimit);

double five_opt_move(instance *inst, int *succ, double curr_cost);

void succ_to_crom(int* succ, int *cromosome, int n);

void crom_to_succ(int *cromosome, int *succ, int cromosome_size, int n);

int update_best(instance *inst, int *succ, double succ_cost);

double get_random_solution(instance *inst, int *succ);