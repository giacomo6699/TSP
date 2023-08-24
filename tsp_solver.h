#ifndef TSP_SOLVER_H_  

#define TSP_SOLVER_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  

#include <cplex.h>  
#include <pthread.h>  

#define VERBOSE				    0		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define EPS_SOL				  1e-5
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ
                                 
//data structures  
typedef struct {   
	
	//input data
	int nnodes; 	
	double *xcoord;
	double *ycoord;

	// parameters 
	int randomseed;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	double best_value;						// best solution available    
	int *best_sol;
	int ncols;
	int integer_costs;
	int verbose;
	int grasp_par;
	int patching_heur;
	double t_start;

} instance;  

/**
 * @brief TSP - greedy solver that uses, as principle for the choice, the least distant nodes starting from a random node
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param grasp_parameter 0 <= x <= 1 number that indicates the probability of choosing the second best choice
 * @param timelimit timelimit of the algorithm
 * @param two_opt 0 -> after the algorithm don't do two_opt; 1 -> do the two_opt after the algorithm found the final solution
 * @return double value containing the cost of the best solution found by the algorithm
 */
double TSPgreedy(instance *inst, int* best_sol, double grasp_parameter, double timelimit, int two_opt);

/**
 * @brief TSP - extramileage solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param covered_nodes_succ array of ints containing the nodes already connected to other nodes in a single cycle
 * @param uncovered_nodes_succ array of ints containing the nodes that were not being added to the cycle yet
 * @param covered_nodes_size number of nodes in covered_nodes_succ
 * @param grasp_parameter 0 <= x <= 1 number that indicates the probability of choosing the second best choice
 * @param timelimit timelimit of the algorithm
 * @param two_opt 0 -> after the algorithm don't do two_opt; 1 -> do the two_opt after the algorithm found the final solution
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double TSPextramileage(instance *inst, int* best_sol, int* covered_nodes_succ, int* uncovered_nodes_succ, int covered_nodes_size, double grasp_parameter, double timelimit, int two_opt);

/**
 * @brief TSP - VNS solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param timelimit timelimit of the algorithm
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double VNS(instance *inst, int *best_sol, double timelimit);

/**
 * @brief TSP - TABU Search solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param tenure_mode costant if 0; step-wise otherwise
 * @param timelimit timelimit of the algorithm
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double TABU_search(instance *inst, int *best_sol, int tenure_mode, double timelimit);

/**
 * @brief TSP - simulated annealing solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param timelimit timelimit of the algorithm
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double simulated_annealing(instance *inst, int *best_sol, double timelimit);

/**
 * @brief TSP - genetic algorithm solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param timelimit timelimit of the algorithm
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double genetic_algorithm(instance *inst, int *best_sol, double timelimit);

/**
 * @brief TSP - benders' loop solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param timelimit timelimit of the algorithm
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double benders_loop(instance *inst, int* best_sol, double timelimit);

/**
 * @brief TSP - "improved" benders' loop solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param timelimit timelimit of the algorithm
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double improved_benders_loop(instance *inst, int *best_sol, double timelimit);

/**
 * @brief TSP - Branch & Cut + Callback solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param timelimit timelimit of the algorithm
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double Branch_and_Cut(instance *inst, int *best_sol, double timelimit);

/**
 * @brief TSP - Hard Fixing or Local Branching solver
 * 
 * @param inst instance of the tsp problem
 * @param best_sol array of ints that will contain the best solution found by this algorithm in "succ" form
 * @param timelimit timelimit of the algorithm
 * @param hf_par parameter to choose which algorithm to use: 0 <= hf_par <= 1 for Hard Fixing and hf_par > 1 for Local Branching
 * @return double value containing the cost of the best solution found by the algorithm 
 */
double branch_and_cut_HF(instance *inst, int *best_sol, double timelimit, double hf_par);

#endif