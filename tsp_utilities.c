#include "tsp_solver.h"
#include <time.h>
#include <cplex.h> 

double random01() { 
	return ((double) rand() / RAND_MAX); 
}

void print_uncovered(int *uncovered, int uncovered_len){
	printf("[");
	for (int i = 0; i < uncovered_len; i++){
		printf("%d, ", uncovered[i]);
	}
	printf("]\n");
}

void print_succ(int *succ, int nnodes){
	printf("[");
	for (int i = 0; i < nnodes; i++){
		printf("%d, ", succ[i]);
	}
	printf("]\n");
}

void print_travel(int *succ, instance *inst){
	FILE *out = fopen("./plot/data.dat", "w");
	int *vis = (int*) malloc(inst->nnodes*sizeof(int));
	for (int i = 0; i < inst->nnodes; i++){
		vis[i] = -1;
	}
	for (int i = 0; i < inst->nnodes; i++){
		if (vis[i] >= 0) continue;
		vis[i] = 1;
		int next = i;
		int done = 0;
		while (!done){
			//printf("\nwriting index %d", next);
			//printf(" that corresponds to %f, %f", inst->xcoord[next], inst->ycoord[next]);
			fprintf(out, "%f %f\r\n", inst->xcoord[next], inst->ycoord[next]);
			next = succ[next];
			vis[next] = 1;
			if (next == i)
				done = 1;
		}
		fprintf(out, "%f %f\r\n", inst->xcoord[next], inst->ycoord[next]);
		fprintf(out, "\r\n", inst->xcoord[next], inst->ycoord[next]);
	}
	//printf("Starting print travel");
    fclose(out);
    system("gnuplot ./plot/commands.txt");
}

void print_instance(instance *inst){
	printf("input_file: %s\n", inst->input_file);
	printf("time limit: %lf\n", inst->timelimit);
	printf("random seed: %d\n", inst->randomseed);
	for (int i = 0; i < inst->nnodes; i++){
		printf("node %4d at coordinates (%8.2f, %8.2f)\n", i, inst->xcoord[i], inst->ycoord[i]);
	}
}


void print_error(const char *err) {
    printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1);
}  

/**
 Depending on the number of nodes we maybe want to have an array with all the distances precomputed instead of computing them each time
 "inline" is useful since the fact of calling the function and passing the parameters can slow down the performances, applying "inline" each time 
 the function is called it will be replaced by the precompiler (as for #define) with the instruction that is within the function (modern compilers do this
 automatically)
inline double cost(int i, int j, instance *inst)
{
	return inst->cst[i*inst->nnodes + j];
}*/

double eucl_dist(int i, int j, instance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	if ( !inst->integer_costs ) return sqrt(dx*dx+dy*dy);
	int dis = sqrt(dx*dx+dy*dy) + 0.499999999; // nearest integer 
	return dis;
}        

double compute_sol_cost(instance *inst, int *succ, int succ_size){
    int *vis = (int*) malloc(succ_size*sizeof(int));
    double cost = 0;
    for (int i = 0; i < succ_size; i++){
		if (succ[i] == -1) vis[i] = 1;
        else vis[i] = -1;
    }
    for (int i = 0; i < succ_size; i++){
        if (vis[i] >= 0) continue;
        int next = i;
        do {
            cost += eucl_dist(next, succ[next], inst);
			vis[next] = 1;
            next = succ[next];
        } while (next != i);
    }
	free(vis);
    return cost;
}

int is_fractional(double x) 						// it works for x in [0,1] only
{
	return ( (x > XSMALL) && (x < 1-XSMALL) );
}    

int is_all_integer(int n, const double *x) 			// it works for x_j in [0,1] only
{
	for ( int j = 0; j < n; j++ ) 
	{
		if ( is_fractional(x[j]) ) return 0; 
	}
	return 1;
}

void double_vector_copy(int n, const double *from, double *to) // vector copy
{
	for ( int j = 0; j < n; j++ ) to[j] = from[j];
}

void choose_starter(int *starter, instance *inst){
	*starter = 0;
}

void choose_random_starter(int *starter, instance *inst){
	*starter = ((double) rand() / RAND_MAX) * inst->nnodes;
}

void choose_starters(int *a, int *b, instance *inst)
{
	double max = 0;
	int max_j = -1;
	int max_i = -1;
	for (int i = 0; i < inst->nnodes - 1; i++){
		for (int j = i+1; j < inst->nnodes; j++){
			double dist = eucl_dist(i, j, inst);
			if (dist > max){
				max = dist;
				max_i = i;
				max_j = j;
			}
		}
	}
	if (max_j == -1 || max_i == -1)
		print_error("Error in computing distances to choose the starters");
	*a = max_i;
	*b = max_j;
}

void choose_random_starters(int *a, int *b, instance *inst)
{
	*a = ((double) rand() / RAND_MAX) * inst->nnodes;
	*b = ((double) rand() / RAND_MAX) * inst->nnodes;
	while (*a == *b){
		*b = ((double) rand() / RAND_MAX) * inst->nnodes;
	}
}

void reverse_edges(int b, int a_prime, int *succ){
	// start from a_prime and change the successor of the following nodes until you reach b
	int first = a_prime;
	int second = succ[a_prime];
	int third = succ[second];
	while (b != first){
		succ[second] = first;
		first = second;
		second = third;
		third = succ[third];
	}
}

double two_opt_move(instance *inst, int *succ, double curr_cost){
	// find the minimum delta_cost
	double min_delta_cost = INFINITY;
	int min_a = -1;
	int min_b = -1;
	for (int a = 0; a < inst->nnodes; a++){
		for (int b = a + 1; b < inst->nnodes; b++){
			//if (a == succ[b] || b == succ[a]) continue; // discard triangles NOT NEEDED, CAN'T HAVE TRIANGLES
			int a_prime = succ[a];
			int b_prime = succ[b];
			double delta_cost = (eucl_dist(a, b, inst) + eucl_dist(a_prime, b_prime, inst)) - (eucl_dist(a, a_prime, inst) + eucl_dist(b, b_prime, inst));
			if (delta_cost < min_delta_cost){
				min_delta_cost = delta_cost;
				min_a = a;
				min_b = b;
			}
		}
	}
	// if the minimum delta_cost is < 0 we have an improvement
	if (min_delta_cost < 0){
		int min_a_prime = succ[min_a];
		int min_b_prime = succ[min_b];
		// invert edges and than remove the cross
		reverse_edges(min_b, min_a_prime, succ);
		succ[min_a] = min_b;
		succ[min_a_prime] = min_b_prime;
		return curr_cost + min_delta_cost;
	} else {
		return curr_cost;
	}
}

double second();

double two_opt_improvement(instance *inst, int *succ, double curr_cost, double timelimit){
	double start = second();
	double end = 0;
	int improvement = 1;
	while(improvement && (end - start) < timelimit){
		double improved_cost = two_opt_move(inst, succ, curr_cost);
		if (improved_cost >= curr_cost - EPS_SOL){
			improvement = 0;
			return curr_cost;
		} else {
			curr_cost = improved_cost;
		}
		end = second();
	}
	return curr_cost;
}

int update_best(instance *inst, int *succ, double succ_cost){
	double computed_cost = compute_sol_cost(inst, succ, inst->nnodes);
	if (fabs(computed_cost - succ_cost) > EPS_SOL) print_error("update_best(): error in solution cost");
	if (succ_cost < inst->best_value - EPS_SOL){
		inst->best_value = succ_cost;
		for (int i = 0; i < inst->nnodes; i++){
			inst->best_sol[i] = succ[i];
		}
		return 1;
	}
	return 0;
}

double five_opt_move(instance *inst, int *succ, double curr_cost){
	int n = inst->nnodes;
	if (n <= 2) print_error("The number of nodes is too small");

	// pick 5 nodes at random such that they are different from each other and not consecutives
	int cand[5];
	int cand_size = 0;
	while(cand_size < 5){
		cand[cand_size] = rand() % n;
		for (int i = 0; i < cand_size; i++){
			if (cand[cand_size] == cand[i] || cand[cand_size] == succ[cand[i]]){
				cand_size--;
				break;
			}
		}
		cand_size++;
	}

	// reorder the candidates following the travel order
	int breaks[5];
	int breaks_succ[5];
	int node = 0;
	int index = 0;
	while(index < 5){
		for (int j = 0; j < 5; j++){
			if (node == cand[j]){
				breaks[index] = node;
				breaks_succ[index] = succ[node];
				index++;
				break;
			}
		}
		node = succ[node];
	}

	// ------ perform the changes and update the cost ------
	double removed_br1 = eucl_dist(breaks[1], breaks_succ[1], inst);
	reverse_edges(breaks[1], breaks_succ[0], succ);

	double removed_br3 = eucl_dist(breaks[3], breaks_succ[3], inst);
	reverse_edges(breaks[3], breaks_succ[2], succ);

	// save the previous distance
	double removed_br0 = eucl_dist(breaks[0], breaks_succ[0], inst);
	// perform the change
	succ[breaks[0]] = breaks[3];
	// save the new distance
	double added_br0 = eucl_dist(breaks[0], succ[breaks[0]], inst);

	succ[breaks_succ[2]] = breaks_succ[3];
	double added_succ_br2 = eucl_dist(breaks_succ[2], succ[breaks_succ[2]], inst);

	double removed_br4 = eucl_dist(breaks[4], breaks_succ[4], inst);
	succ[breaks[4]] = breaks[1];
	double added_br4 = eucl_dist(breaks[4], succ[breaks[4]], inst);

	succ[breaks_succ[0]] = breaks_succ[1];
	double added_succ_br0 = eucl_dist(breaks_succ[0], succ[breaks_succ[0]], inst);

	double removed_br2 = eucl_dist(breaks[2], breaks_succ[2], inst);
	succ[breaks[2]] = breaks_succ[4];
	double added_br2 = eucl_dist(breaks[2], succ[breaks[2]], inst);

	double new_edges_cost = added_br0 + added_br2 + added_br4 + added_succ_br0 + added_succ_br2;
	double old_edges_cost = removed_br0 + removed_br1 + removed_br2 + removed_br3 + removed_br4;

	curr_cost = curr_cost + new_edges_cost - old_edges_cost;
	// ------ perform the changes and update the cost ------

	return curr_cost;
}

void succ_to_crom(int* succ, int *cromosome, int n){
	int next = 0;
	for (int i = 0; i < n; i++){
		cromosome[i] = succ[next];
		next = succ[next];
	}
}

void crom_to_succ(int *cromosome, int *succ, int cromosome_size, int n){
	for (int i = 0; i < n; i++){
		succ[i] = -1;
	}
	for (int i = 0; i < cromosome_size - 1; i++){
		succ[cromosome[i]] = cromosome[i + 1];
	}
	succ[cromosome[cromosome_size-1]] = cromosome[0];
}

double get_random_solution(instance *inst, int *succ){
	int n = inst->nnodes;
	double cost = 0;

	int *uncovered = (int*) calloc(n - 1, sizeof(int));
	int uncovered_length = n - 1;
	for (int i = 0; i < uncovered_length; i++){
		uncovered[i] = i + 1;
	}

	int curr_covered = 0;
	while (uncovered_length > 0) {
		int index = rand() % uncovered_length;
		succ[curr_covered] = uncovered[index];
		cost += eucl_dist(curr_covered, uncovered[index], inst);
		curr_covered = uncovered[index];
		uncovered[index] = uncovered[--uncovered_length];
	}
	cost += eucl_dist(0, curr_covered, inst);
	succ[curr_covered] = 0;

	free(uncovered);
	return cost;
}