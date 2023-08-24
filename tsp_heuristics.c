#include "tsp_solver.h"
#include "tsp_utilities.h"
#include <time.h>
#include <cplex.h> 
#include "tsp_utilities.h"
#include "chrono.h"

// ---------------------------- HEURISTICS --------------------------------

double TSPgreedy(instance *inst, int* best_sol, double grasp_parameter, double timelimit, int two_opt)
{
	double start = second();
	double end = 0;

	double greedy_timelimit = timelimit;
	double two_opt_timelimit = 0;

	if (two_opt){
		greedy_timelimit = timelimit/5*4;
		two_opt_timelimit = timelimit/5;
	}

	int *succ = (int*) malloc(inst->nnodes*sizeof(int));
	int *uncovered = (int*) malloc((inst->nnodes - 1) * sizeof(int));
	double best_cost = INFINITY;

	while (end - start < greedy_timelimit){
		// get a starter
		int starter;
		choose_random_starter(&starter, inst);

		// i is the covered node that we are considering
		int i = starter;

		// uncovered_len is nnodes - 1 since the starter is already taken out of uncovered
		int uncovered_len = (inst->nnodes - 1);
		// sol_cost contains the sum of the distances in the travel
		double sol_cost = 0;

		// initialize the uncovered vector
		for (int h = 0; h < uncovered_len; h++){
			if (h >= starter){
				uncovered[h] = h + 1;
			} else {
				uncovered[h] = h;
			}
		}

		// while there are still some nodes that are uncovered
		while(uncovered_len > 0){
			// min_dist will contain the minimum distance between i and the uncovered nodes
			double min_dist = INFINITY;
			// second_min_dist will contain the second minimum distance between i and the uncovered nodes
			double second_min_dist = INFINITY;
			// min_j will contain the uncovered node that is at distance min_dist from i
			int min_j = -1;
			// second_min_j will contain the uncovered node that is at distance second_min_dist from i
			int second_min_j = -1;
			// min_j_index will contain the index of min_j inside of the uncovered vector (needed to do the resize operation)
			int min_j_index = -1;
			// second_min_j_index will contain the index of second_min_j inside of the uncovered vector (needed to do the resize operation)
			int second_min_j_index = -1;

			// random value
			double grasp_value = random01();
			//printf("GRASP VALUE %f\n", grasp_value);

			// iterate through all the uncovered nodes
			for (int index = 0; index < uncovered_len; index++){
				int j = uncovered[index];
				// compute the distance between i and the uncovered node
				double curr_dist = eucl_dist(i, j, inst);
				// if curr_dist is the new minimum update the variables
				if (curr_dist < min_dist){
					min_dist = curr_dist;
					min_j = j;
					min_j_index = index;
				} else if (curr_dist < second_min_dist && grasp_value < grasp_parameter){ // if curr_dist is not the new minimum but it's the new second minimum update the related variables
					second_min_dist = curr_dist;
					second_min_j = j;
					second_min_j_index = index;
				} 
			}
			// if for checking that at least one time the variables have been update
			if (min_j == -1 || min_j_index == -1)
				print_error("ERROR IN TSPgreedy");

			// check the grasp_value --> randomizing the choice taking as successor the second best choice
			if (grasp_value < grasp_parameter && second_min_j != -1){ // check the probability and also that there were at least 2 nodes among the uncovered ones
				// second best choice
				// remove second_min_j from the uncovered vector
				uncovered[second_min_j_index] = uncovered[--uncovered_len];
				// add second_min_j as i's successor
				succ[i] = second_min_j;
				// update the value of i since we have to consider the next covered node
				i = second_min_j;
				// update the cost of the solution
				sol_cost += second_min_dist;
				if (uncovered_len == 0){
					succ[second_min_j] = starter;
					sol_cost += eucl_dist(second_min_j, starter, inst);
				}
			} else {
				// first best choice
				// remove min_j from the uncovered vector
				uncovered[min_j_index] = uncovered[--uncovered_len];
				// add min_j as i's successor
				succ[i] = min_j;
				// update the value of i since we have to consider the next covered node
				i = min_j;
				// update the cost of the solution
				sol_cost += min_dist;
				if (uncovered_len == 0){
					succ[min_j] = starter;
					sol_cost += eucl_dist(min_j, starter, inst);
				}
			}
		}
		
		if (sol_cost < best_cost - EPS_SOL){
			best_cost = sol_cost;
			for (int i = 0; i < inst->nnodes; i++){
				best_sol[i] = succ[i];
			}
		}

		end = second();
	}


	if (two_opt)
		best_cost = two_opt_improvement(inst, best_sol, best_cost, two_opt_timelimit);
	update_best(inst, best_sol, best_cost);
	if (inst->verbose >= 60)
		printf("\nGreedy Solution found in %f seconds has cost: %f\n", second() - start, best_cost);
	free(uncovered);
	free(succ);
	return best_cost;
}

double TSPextramileage(instance *inst, int* best_sol, int* covered_nodes_succ, int* uncovered_nodes_succ, int covered_nodes_size, double grasp_parameter, double timelimit, int two_opt)
{
	double start = second();
	double end = 0;

	double extram_timelimit = timelimit;
	double two_opt_timelimit = 0;
	if (two_opt){
		extram_timelimit = timelimit/5*4;
		two_opt_timelimit = timelimit/5;
	}


	double best_cost = INFINITY;
	int uncovered_len = inst->nnodes - covered_nodes_size;
	double sol_cost = 0;

	while (end - start < extram_timelimit){
		
		int a;

		if (covered_nodes_size == 0){

			// get the starters
			int b;
			choose_starters(&a, &b, inst);
			// uncovered_len is nnodes - 2 since the starters are already taken out of uncovered
			uncovered_len = inst->nnodes - 2;
			// sol_cost contains the sum of the distances in the travel
			sol_cost = 2 * eucl_dist(a, b, inst);

			// set the starting point by adding the starters successors
			covered_nodes_succ[a] = b;
			covered_nodes_succ[b] = a;

			// initialize the uncovered vector skipping the starters (that are covered)
			for (int h = 0; h < uncovered_len; h++)
			{
				if (h >= (a - 1) && h >= (b - 1))
				{
					uncovered_nodes_succ[h] = h + 2;
				}
				else if (h >= a && h < b)
				{
					uncovered_nodes_succ[h] = h + 1;
				}
				else if (h < a && h >= b)
				{
					uncovered_nodes_succ[h] = h + 1;
				}
				else
				{
					uncovered_nodes_succ[h] = h;
				}
			}
		} else {
			for (int i = 0; i < inst->nnodes; i++){
				if (covered_nodes_succ[i] != -1){
					a = covered_nodes_succ[i];
					break;
				}
			}
			sol_cost = compute_sol_cost(inst, covered_nodes_succ, inst->nnodes);
		}

		while(uncovered_len != 0){ // while there are still some nodes that are uncovered
			double min_cost = INFINITY; // will contain the minimum cost found for the pair (i, h)
			int min_h = -1; // will contain the uncovered node that has min_cost with i
			int min_i = -1; // will contain the uncovered node that has min_cost with h
			int uncovered_index_h = -1; // will contain the index of min_h inside of the uncovered vector (needed to do the resize operation)
			double second_min_cost = INFINITY; // will contain the second minimum cost found for the pair (i, h)
			int second_min_h = -1; // will contain the uncovered node that has second_min_cost with i
			int second_min_i = -1; // will contain the uncovered node that has second_min_cost with h
			int second_uncovered_index_h = -1; // will contain the index of second_min_h inside of the uncovered vector (needed to do the resize operation)
			int i = -1; // is the considered node (among the covered ones)

			double grasp_value = random01(); // random value

			// while i != a is for detecting when we finished considering the covered nodes (we ended up at the end of the travel)
			// so we will iterate through all the covered nodes without using an additional structure but only by using the succ vector
			while (i != a){
				if (i == -1) i = a; // if we are at the starting iteration put i equal to the starter (done here cause of the while condition)
				for (int h = 0; h < uncovered_len; h++){ // iterate through all the uncovered nodes
					// compute the delta between i and uncovered[h]
					double delta = eucl_dist(i, uncovered_nodes_succ[h], inst) + eucl_dist(uncovered_nodes_succ[h], covered_nodes_succ[i], inst) - eucl_dist(i, covered_nodes_succ[i], inst);
					// if delta is the new minimum, update the variables
					if (delta < min_cost){
						min_cost = delta;
						min_h = uncovered_nodes_succ[h];
						uncovered_index_h = h;
						min_i = i;
					} else if (delta < second_min_cost && grasp_value < grasp_parameter){ // if delta is not the new minimum but it's the new second minimum update the related variables
						second_min_cost = delta;
						second_min_h = uncovered_nodes_succ[h];
						second_uncovered_index_h = h;
						second_min_i = i;
					}
				}
				// update i to be the next node among the covered ones (following the successors order)
				i = covered_nodes_succ[i];
			}
			// check if at least one time min_cost have been updated
			if (min_h == -1)
				print_error("Error in computing the minimum delta");
			
			// check the grasp_value --> randomizing the choice taking as successor the second best choice
			if (grasp_value < grasp_parameter && second_min_h != -1){
				covered_nodes_succ[second_min_h] = covered_nodes_succ[second_min_i]; // update the successor of second_min_h to be the old successor of second_min_i
				covered_nodes_succ[second_min_i] = second_min_h; // update the successor of second_min_i to be second_min_h
				uncovered_nodes_succ[second_uncovered_index_h] = uncovered_nodes_succ[--uncovered_len]; // update the uncovered vector
				sol_cost += second_min_cost; // update the cost of the solution
			} else {
				covered_nodes_succ[min_h] = covered_nodes_succ[min_i]; // update the successor of min_h to be the old successor of min_i
				covered_nodes_succ[min_i] = min_h; // update the successor of min_i to be min_h
				uncovered_nodes_succ[uncovered_index_h] = uncovered_nodes_succ[--uncovered_len]; // update the uncovered vector
				sol_cost += min_cost; // update the cost of the solution
			}
		}

		if (sol_cost < best_cost - EPS_SOL){
			best_cost = sol_cost;
			for (int i = 0; i < inst-> nnodes; i++){
				best_sol[i] = covered_nodes_succ[i];
			}
		}

		if (covered_nodes_size != 0){ // only one iteration, no need for grasp in case of a partial solution as input
			break;
		} else {
			end = second();
		}
	}
	if (two_opt)
		best_cost = two_opt_improvement(inst, best_sol, best_cost, two_opt_timelimit);
	update_best(inst, best_sol, best_cost);
	if (inst->verbose >= 60)
		printf("\nExtra Mileage Solution found in %f seconds has cost: %f\n", second() - start, best_cost);
	return best_cost;
}


// ---------------------------- METAHEURISTICS --------------------------------


double VNS(instance *inst, int *best_sol, double timelimit){
	int n = inst->nnodes;
	double start = inst->t_start;
	double end = 0;

	double greedy_timelimit = timelimit/5;
	
	int *succ = (int*) malloc(n * sizeof(int));

	// call a solving method + 2-opt
	double curr_cost = TSPgreedy(inst, succ, 0.1, greedy_timelimit, 1);
	double best_cost = curr_cost;
	for (int i = 0; i < n; i++){
		best_sol[i] = succ[i];
	}

	//FILE *out = fopen("./plot/costs.txt", "w");
	//fprintf(out, "$$%f\n", curr_cost);

	while (end - start < timelimit){
		// kick
		curr_cost = five_opt_move(inst, succ, curr_cost);

    	//fprintf(out, "$$%f\n", curr_cost);

		double remaining_tl = timelimit - (second() - start);

		curr_cost = two_opt_improvement(inst, succ, curr_cost, remaining_tl);

    	//fprintf(out, "$$%f\n", curr_cost);

		if (curr_cost < best_cost - EPS_SOL){
			best_cost = curr_cost;
			for (int i = 0; i < n; i++){
				best_sol[i] = succ[i];
			}
		}

		end = second();
	}

	//fclose(out);

	update_best(inst, best_sol, best_cost);

	if (inst->verbose >= 10)
		printf("\nVNS Search Solution found in %f seconds has cost: %f\n", second() - start, best_cost);

	free(succ);
	
	return best_cost;

}

double TABU_search(instance *inst, int *best_sol, int tenure_mode, double timelimit){
	double start = inst->t_start;
	double end = 0;
	int n = inst->nnodes;

	double greedy_timelimit = timelimit/5;

	int *succ = (int*) malloc(n * sizeof(int));
	int *tabu_list = (int*) malloc(n * sizeof(int));
	for(int i = 0; i < n; i++){
		tabu_list[i] = -999999999;
	}

	// ------ call a solving method + 2-opt ------
	double curr_cost = TSPgreedy(inst, succ, 0.1, greedy_timelimit, 1);
	double best_cost = curr_cost;
	for (int i = 0; i < n; i++){
		best_sol[i] = succ[i];
	}

	int iter_index = 0;
	while(end - start < timelimit){
		int tenure;
		if (tenure_mode == 0){
			tenure = 100;
		} else {
			tenure = iter_index / 100 % 2 ? 20 : 100;
		}
		double min_delta_cost = INFINITY;
		int min_a = -1;
		int min_b = -1;

		// ------ find the best move not involving the tabu nodes ------
		for (int a = 0; a < n; a++){
			for (int b = a + 1; b < n; b++){
				if (a == succ[b] || b == succ[a]) continue; // discard triangles
				int a_prime = succ[a];
				int b_prime = succ[b];
				double delta_cost = (eucl_dist(a, b, inst) + eucl_dist(a_prime, b_prime, inst)) - (eucl_dist(a, a_prime, inst) + eucl_dist(b, b_prime, inst));
				if (iter_index - tabu_list[a] > tenure && iter_index - tabu_list[b] > tenure){
					if (delta_cost < min_delta_cost){
						min_delta_cost = delta_cost;
						min_a = a;
						min_b = b;
					}
				} else {
					if (curr_cost + delta_cost < best_cost){ // check for aspiration criteria
						min_delta_cost = delta_cost;
						min_a = a;
						min_b = b;
					}
				}
			}
		}

		// ------ perform the changes in the graph ------
		int min_a_prime = succ[min_a];
		int min_b_prime = succ[min_b];
		reverse_edges(min_b, min_a_prime, succ);
		succ[min_a] = min_b;
		succ[min_a_prime] = min_b_prime;

		// update the cost
		curr_cost += min_delta_cost;

		// set tabu list values
		tabu_list[min_a] = iter_index;
		tabu_list[min_b] = iter_index;

		// update the best cost and the best solution
		if (curr_cost < best_cost - EPS_SOL){
			best_cost = curr_cost;
			for (int i = 0; i < n; i++){
				best_sol[i] = succ[i];
			}
		}

		iter_index++;
		end = second();
	}

	update_best(inst, best_sol, best_cost);

	if (inst->verbose >= 30)
		printf("\nTabu Search Solution found in %f seconds has cost: %f\n", second() - start, best_cost);

	free(tabu_list);
	free(succ);

	return best_cost;
}

// WORKS BUT NOT SURE OF THE PARAMETERS DEFINITION
double simulated_annealing(instance *inst, int *best_sol, double timelimit){
	int n = inst->nnodes;
	double t_start = inst->t_start;
	double t_end = 0;

	int *succ = (int*) malloc(n * sizeof(int));
	double curr_cost = TSPgreedy(inst, succ, 0.0, 0.0, 0); // greedy without two_opt
	double scaler_cost = curr_cost; // for my implementation can be removed
	
	double best_cost = curr_cost;
	for (int i = 0; i < n; i++){
		best_sol[i] = succ[i];
	}

	double old_min = 0;
	double old_max = timelimit;
	while(t_end - t_start < timelimit){
		// compute delta_cost of hypotethical random two-opt
		int a = rand() % n;
		int b = rand() % n;
		while (a == b){
			b = rand() % n;
		}
		int a_prime = succ[a];
		int b_prime = succ[b];
		double delta_cost = (eucl_dist(a, b, inst) + eucl_dist(a_prime, b_prime, inst)) - (eucl_dist(a, a_prime, inst) + eucl_dist(b, b_prime, inst));
		
		if (delta_cost <= 0){
			// accepted: make the changes to succ
			reverse_edges(b, a_prime, succ);
			succ[a] = b;
			succ[a_prime] = b_prime;
			curr_cost += delta_cost;
		} else {
			double old_value = timelimit - (second() - t_start);
			double new_max = - delta_cost / (scaler_cost * log(1/2)); // in order to have 50% probability
			double new_min = delta_cost / (scaler_cost * 1000); // in order to have almost 0% probability
			double temperature = ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min; // transform from [0, timelimit] range to [new_min, new_max] range
			double accept_prob = exp(-delta_cost/(scaler_cost * temperature));
			if (random01() < accept_prob){
				// accepted: make the changes to succ
				reverse_edges(b, a_prime, succ);
				succ[a] = b;
				succ[a_prime] = b_prime;
				curr_cost += delta_cost;
			}
		}

		if (curr_cost < best_cost - EPS_SOL){
			best_cost = curr_cost;
			for (int i = 0; i < n; i++){
				best_sol[i] = succ[i];
			}
		}

		t_end = second();
	}
	update_best(inst, best_sol, best_cost);
	free(succ);
	return best_cost;
}

double apply_mutation(instance *inst, int *cromosome, int cromosome_size, int mutations_number, double curr_cost){
	int changes_number = rand() % mutations_number + 1; // at least one mutation, max mutations_nummber mutations
	for (int i = 0; i < changes_number; i++){
		int first_index = rand() % cromosome_size;
		int second_index = rand() % cromosome_size;
		while (first_index == second_index)
			second_index = rand() % cromosome_size;

		int first_previous_index = (first_index - 1 >= 0) ? first_index - 1 : cromosome_size - 1;
		int first_next_index = (first_index + 1 < cromosome_size) ? first_index + 1 : 0;
		double cost_to_remove = eucl_dist(cromosome[first_previous_index], cromosome[first_index], inst) + eucl_dist(cromosome[first_index], cromosome[first_next_index], inst);

		int second_previous_index = (second_index - 1 >= 0) ? second_index - 1 : cromosome_size - 1;
		int second_next_index = (second_index + 1 < cromosome_size) ? second_index + 1 : 0;
		cost_to_remove += eucl_dist(cromosome[second_previous_index], cromosome[second_index], inst) + eucl_dist(cromosome[second_index], cromosome[second_next_index], inst);

		double cost_to_add = eucl_dist(cromosome[first_previous_index], cromosome[second_index], inst) + eucl_dist(cromosome[second_index], cromosome[first_next_index], inst);
		cost_to_add += eucl_dist(cromosome[second_previous_index], cromosome[first_index], inst) + eucl_dist(cromosome[first_index], cromosome[second_next_index], inst);

		int temp = cromosome[first_index];
		cromosome[first_index] = cromosome[second_index];
		cromosome[second_index] = temp;

		curr_cost = curr_cost + cost_to_add + cost_to_remove;
	}

	return curr_cost;
}

// NOT COMPLETELY DEBUGGED, COULD CONTAIN BUGS 
double genetic_algorithm(instance *inst, int *best_sol, double timelimit)
{
	double start = inst->t_start;
	double end = 0;
	double generation_tl = timelimit/5;

	int n = inst->nnodes;
	int pop_size = 1000;
	int max_size = pop_size + pop_size / 4 + pop_size / 10;
	int *population = (int *)malloc(max_size * n * sizeof(int));
	double *population_costs = (double *) malloc(max_size * sizeof(double));
	int *succ = (int *)malloc(n * sizeof(int));
	double best_cost = INFINITY;

	FILE *out = fopen("./plot/costs.txt", "w");


	// ------ generation of the population ------ popsize elements + popsize/4 children + popsize/10 mutated
	// generate 90% of population using grasp ad 10% in a random way then apply two_opt to all
	for (int i = 0; i < 1000; i++)
	{
		// used this remaining_tl because inside TSPgreedy a fraction of the timelimit is dedicated to grasp and another fraction to two_opt
		// since two_opt may not use all the timelimit reserved to it, we don't want to waste that time and so we reallocate it among all the remaining future instances
		double remaining_tl = generation_tl - (second() - start);
		if (i < 900)
			population_costs[i] = TSPgreedy(inst, succ, 0.2, remaining_tl/(1000 - i), 1);
		else {
			double random_cost = get_random_solution(inst, succ);
			random_cost = two_opt_improvement(inst, succ, random_cost, remaining_tl/(1000 - i));
			population_costs[i] = random_cost;
			update_best(inst, succ, random_cost);
		}
		succ_to_crom(succ, &population[i * n], n); // convert cromosome and paste it into the population
	}
	
	if (inst->verbose >= 50) printf("\nInitial population generated\n");

	int *nodes_flags = (int*) malloc(n * sizeof(int));
	int iterations = 0;

	while (end - start < timelimit)
	{
		iterations++;

		// ------ crossover and repairment ------
		for (int i = 0; i < pop_size / 4; i++) // for popsize/4 times take 2 elements at random (the parents)
		{
			int parent1_index = rand() % pop_size; // random number in [0, popsize - 1]
			int parent2_index = rand() % pop_size;
			while (parent1_index == parent2_index)
				parent2_index = rand() % pop_size;
			int cutting_point = rand() % (3 * n / 4 - n / 4 + 1) + n / 4; // for each pair set the cutting point picking a random value between n/4 and 3n/4
			// create a child appending to the population array the first part of parent1 and the second one of parent2
			// we skip the nodes that have been already added and add the missing nodes to an uncover vector
			int child_size = 0;
			for (int j = 0; j < n; j++){
				nodes_flags[j] = 0;
			}
			for (int j = 0; j < n; j++)
			{
				int child_elem_index = (pop_size + i) * n + child_size;
				int parent_elem_index;
				if (j < cutting_point)
				{
					parent_elem_index = parent1_index * n + j;
				}
				else
				{
					parent_elem_index = parent2_index * n + j;
				}
				if (nodes_flags[population[parent_elem_index]] == 0)
				{
					population[child_elem_index] = population[parent_elem_index];
					nodes_flags[population[parent_elem_index]] = 1;
					child_size++;
				}
			}
			// collect the missing nodes inside the array uncovered_nodes
			int *uncovered_nodes = (int*) malloc((n - child_size) * sizeof(int));
			int uncovered_index = 0;
			for (int j = 0; j < n; j++)
			{
				if (nodes_flags[j] == 0)
				{
					uncovered_nodes[uncovered_index] = j;
					uncovered_index++;
				}
			}

			int child_index = (pop_size + i) * n;
			crom_to_succ(&population[child_index], succ, child_size, n); // convert from the cromosome representation to the successors representation

			int *best_extram = (int*) malloc(n * sizeof(int));
			population_costs[pop_size + i] = TSPextramileage(inst, best_extram, succ, uncovered_nodes, child_size, 0.0, 0.0, 0); // pass it to extra_mileage
			population_costs[pop_size + i] = two_opt_improvement(inst, best_extram, population_costs[pop_size + i], INFINITY); // full 2-opt
			succ_to_crom(best_extram, &population[child_index], n); // reconvert it to add the final child to the population
			free(best_extram);
			free(uncovered_nodes);
		}
		if (inst->verbose >= 60) printf("\nOffspring added to the population\n");

		// ------ mutation ------
		int *mutated_elems = (int*) malloc(pop_size / 10 * sizeof(int));
		for (int i = 0; i < pop_size / 10; i++)
		{
			int mutated_index;
			int already_present;
			do
			{
				mutated_index = rand() % pop_size;
				already_present = 0;
				for (int j = 0; j < i; j++)
				{
					if (mutated_elems[j] == mutated_index)
					{
						already_present = 1;
					}
				}

			} while (already_present);
			mutated_elems[i] = mutated_index;
			int mutated_child_index = (pop_size + pop_size / 4 + i) * n;
			for (int j = 0; j < n; j++)
			{
				population[mutated_child_index + j] = population[mutated_index * n + j];
			}
			population_costs[(pop_size + pop_size / 4 + i)] = apply_mutation(inst, &population[mutated_child_index], n, 3, population_costs[mutated_index]);
		}
		free(mutated_elems);
		if (inst->verbose >= 60) printf("\nMutation phase completed\n");

		// ------ find the champion ------
		double min_cost = INFINITY;
		double max_cost = -INFINITY;
		int champion_index = -1;
		for (int i = 0; i < max_size; i++)
		{
			if (population_costs[i] < min_cost)
			{
				min_cost = population_costs[i];
				champion_index = i;
			}
			if (population_costs[i] > max_cost)
			{
				max_cost = population_costs[i];
			}
		}
		if (min_cost < best_cost){
			best_cost = min_cost;
			crom_to_succ(&population[champion_index * n], best_sol, n, n);
		}
		fprintf(out, "$$%f\n", min_cost);
		update_best(inst, best_sol, best_cost);
		if (inst->verbose >= 50) printf("\nCurrent champion has index %d and cost %f\n", champion_index, min_cost);

		// ------ kill popsize/4 + popsize/10 individuals ------
		int kill_count = 0;
		int new_pop_size = max_size;
		while (kill_count < pop_size / 4 + pop_size / 10){
			int index = rand() % new_pop_size; // take a number at random in [0, new_pop_size - 1]
			if (index == champion_index)
				continue;
			double old_value = population_costs[index];
			double prob = ((old_value - min_cost) / (max_cost - min_cost)); // normalize the cost of that individual
			if (random01() < prob){	// take another random number and if it is < than the normalized cost of the individual
				kill_count++; // increase the kill_count
				population_costs[index] = population_costs[--new_pop_size];
				population[index * n] = population[new_pop_size * n]; // remove the individual from the population and from the population_costs
			}
		}
		if (inst->verbose >= 50) printf("\nKilling phase completed\n");

		end = second();
	}
	fclose(out);
	if (inst->verbose >= 20) printf("\nNumber of iterations: %d\n", iterations);
	update_best(inst, best_sol, best_cost);
	if (inst->verbose >= 10) printf("\nGenetic algorithm found a solution of cost %f in %f seconds\n", best_cost, end - start);
	free(nodes_flags);
	free(succ);
	free(population_costs);
	free(population);
	return best_cost;
}