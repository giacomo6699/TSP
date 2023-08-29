#include "tsp_solver.h"
#include "tsp_utilities.h"
#include <time.h>
#include <cplex.h> 
#include "chrono.h"

int xpos(int i, int j, instance *inst)                                         
{ 
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}

void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
{    

	double zero = 0.0;  
	char binary = 'B'; 

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

	// add binary var.s x(i,j) for i < j  

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = eucl_dist(i,j,inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s");
		}
	} 

	// add the degree constraints 

	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes, sizeof(double));

	for ( int h = 0; h < inst->nnodes; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1);   
		int nnz = 0;
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	if ( VERBOSE >= 100 ) CPXwriteprob(env, lp, "model.lp", NULL);   
}

//#define DEBUG    // comment out to avoid debugging 
#define EPS 1e-5

void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp)
{   

#ifdef DEBUG
	int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			int k = xpos(i,j,inst);
			if ( fabs(xstar[k]) > EPS && fabs(xstar[k]-1.0) > EPS ) print_error(" wrong xstar in build_sol()");
			if ( xstar[k] > 0.5 ) 
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for (int i = 0; i < inst->nnodes; i++)
	{
		if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
	}	
	free(degree);
#endif
	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		succ[i] = -1;
		comp[i] = -1;
	}
	
	for ( int start = 0; start < inst->nnodes; start++ )
	{
		if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done )  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;  // last arc to close the cycle
		
		// go to the next component...
	}
}

void add_one_sec(instance *inst, CPXENVptr env, CPXLPptr lp, int comp_index, int *comp, CPXCALLBACKCONTEXTptr context){
	int ncols = inst->ncols;
	int *index = (int*) calloc(ncols, sizeof(int));
	double *coeff = (double*) calloc(ncols, sizeof(double));
	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));

	int nnz = 0;
	double rhs = -1;
	char sense = 'L';
	sprintf(cname[0], "sec(%d)", comp_index);
	for (int i = 0; i < inst->nnodes; i++){
		if (comp[i] != comp_index) continue;
		rhs++;
		for (int j = i + 1; j < inst->nnodes; j++){
			if (comp[j] != comp_index) continue;
			index[nnz] = xpos(i, j, inst);
			coeff[nnz] = 1.0;
			nnz++;
		}
	}
	int izero = 0;
	if (context == NULL && (env == NULL || lp == NULL)){
		print_error("Error in add_one_sec arguments");
	} else {	
		if (context == NULL){
			if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, coeff, NULL, &cname[0]))
				print_error("CPXaddrows(): error 1");
		} else {
			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, coeff))
				print_error("CPXcallbackrejectcandidate() error"); // reject the solution and adds one cut 
		}

	}
	
	free(cname[0]);
	free(cname);
	free(coeff);
	free(index);
}

void add_invalid_cons(instance *inst, CPXENVptr env, CPXLPptr lp, double *xh, int K){
	int ncols = inst->ncols;
	int *index = (int*) calloc(ncols, sizeof(int));
	double *coeff = (double*) calloc(ncols, sizeof(double));
	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));

	int nnz = 0;
	double rhs = inst->nnodes - K;
	char sense = 'G';
	sprintf(cname[0], "Local Branch");
	for (int i = 0; i < inst->ncols; i++){
		if (xh[i] > 0.5){
			index[nnz] = i;
			coeff[nnz] = 1.0;
			nnz++;
		}
	}
	int izero = 0;
	if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, coeff, NULL, &cname[0]))
		print_error("CPXaddrows(): error 1");

	free(cname[0]);
	free(cname);
	free(coeff);
	free(index);
}

double benders_loop(instance *inst, int* best_sol, double timelimit)
{  

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);
	
	// Cplex's parameter setting
	if (inst->verbose >= 50) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	else CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);	
	CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit); 

	error = CPXmipopt(env,lp);
	if ( error ) 
	{
		printf("CPX error code %d\n", error);
		print_error("CPXmipopt() error"); 
	}
	
	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;
	double *xstar = (double *) calloc(ncols, sizeof(double));


	if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) print_error("first CPXgetx() error");	


	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int ncomp = 0;
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));

	build_sol(xstar, inst, succ, comp, &ncomp);

	// if there are multiple connected components
	while (ncomp > 1 && timelimit > (second() - inst->t_start)){
		// add sec constraint
		for (int k = 1; k <= ncomp; k++){
			add_one_sec(inst, env, lp, k, comp, NULL);
		}
		CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit - (second() - inst->t_start)); 
		// solve current model
		if ( CPXmipopt(env,lp) ) print_error("CPXmipopt() error after sec"); 
		// retrieve x* solution of current model
		if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) printf("after sec CPXgetx() error");
		// detect the different connected components building the array of the different components and succ
		build_sol(xstar, inst, succ, comp, &ncomp);
	}

	for (int i = 0; i < inst->nnodes; i++){
		best_sol[i] = succ[i];
	}
	double curr_cost = 0;
	CPXgetobjval(env, lp, &curr_cost);

	free(comp);
	free(xstar);
	free(succ);
	
	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	return curr_cost;

}

double patching_heuristic(instance *inst, CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int *ncomp, double curr_cost, CPXCALLBACKCONTEXTptr context){
	double initial_cost = curr_cost;
	int* temp = (int*) malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++){
		temp[i] = succ[i];
	}
	while(*ncomp > 1 && (second() - inst->t_start) < inst->timelimit){
		double min_cost = INFINITY;
		int a = -1;
		int b = -1;
		for (int i = 0; i < inst->nnodes; i++){
			for (int j = 0; j < inst->nnodes; j++){
				if (comp[i] < comp[j]){
					double old_cost = eucl_dist(i, temp[i], inst) + eucl_dist(j, temp[j], inst);
					double new_cost = eucl_dist(i, temp[j], inst) + eucl_dist(j, temp[i], inst);
					if ((new_cost - old_cost) < min_cost) {
						min_cost = new_cost - old_cost;
						a = i;
						b = j;
					}
				}
			}
		}
		int node = b;
		int compA = comp[a];
		do {
			comp[node] = compA;
			node = temp[node];
		}
		while (node != b);
		
		(*ncomp)--;

		int a_prime = temp[a];
		temp[a] = temp[b];
		temp[b] = a_prime;

		curr_cost += min_cost;

		if (*ncomp > 1) {
			if (context == NULL && (env == NULL || lp == NULL))
				print_error("Error in patching_heuristic arguments");
			else {	
				if (context == NULL){
					add_one_sec(inst, env, lp, compA, comp, NULL);
				} else {
					add_one_sec(inst, NULL, NULL, compA, comp, context);
				}
			}
		}
	}
	if (*ncomp > 1){
		free(temp);
		return -1;
	} else {
		for (int i = 0; i < inst->nnodes; i++){
			succ[i] = temp[i];
		}
		free(temp);
		double cost = two_opt_improvement(inst, succ, curr_cost, inst->timelimit - (second() - inst->t_start));
		return cost;
	}
}

double improved_benders_loop(instance *inst, int *best_sol, double timelimit){
	double t_start = inst->t_start;
	double t_end = 0;
	
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);
	
	if (inst->verbose >= 50) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	else CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);	

	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;
	double *xstar = (double *) calloc(ncols, sizeof(double));

	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int ncomp = 0;
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));

	double best_obj_val = -9999999;
	double LB = -9999999;
	double UB = INFINITY;
	int iteration = 0;

	// if there are multiple connected components
	while (LB < 0.9999 * UB && (t_end - inst->t_start) < timelimit){
		CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit - (second() - t_start));
		CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);
		
		error = CPXmipopt(env,lp); // solve the problem
		if ( error ) 
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}

		if (CPXgetbestobjval(env, lp, &best_obj_val)) print_error("CPXgetbestobjval() error"); // get the best objective value

		double curr_cost;
		CPXgetobjval(env, lp, &curr_cost);

		LB = (best_obj_val > LB) ? best_obj_val : LB; // if the best objective value is greater than the known bound, update it

		if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) printf("CPXgetx() did not find a better solution");	 // get the solution

		build_sol(xstar, inst, succ, comp, &ncomp); // convert from xstar to succ and get the components informations

		if (ncomp > 1){ // if there are more than one connected component -> add the correspondent constraints
			for (int k = 1; k <= ncomp; k++){ add_one_sec(inst, env, lp, k, comp, NULL); }
		}

		if (inst->verbose >= 60)
			printf("\nCost before patching: %f", curr_cost);
		
		double patch_result = patching_heuristic(inst, env, lp, succ, comp, &ncomp, curr_cost, NULL); // apply patching heuristic + two_opt
		if (patch_result < 0){
			printf("\nPatching Heuristic exceeded the time limit");
		} else {
			curr_cost = patch_result; // in this case the function patching_heuristic has changed the values in the succ vector
		}
		if (inst->verbose >= 60)
			printf("\nCost after patching: %f\n", curr_cost);

		if (curr_cost > 0 && curr_cost < UB - EPS_SOL){ // if we found a better solution, update the upper bound and the best_sol
			UB = curr_cost;
			for (int i = 0; i < inst->nnodes; i++){
				best_sol[i] = succ[i];
			}
		}

		if (inst->verbose >= 30)
			printf("Iteration %d: lb = %11.3f, incumbent = %11.3f, total time = %7.2f\n", iteration, LB, UB, second() - inst->t_start);
		iteration++;
		t_end = second();
	}

	update_best(inst, best_sol, UB);

	free(comp);
	free(succ);
	free(xstar);
	
	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	return UB;
}



static int CPXPUBLIC incumbent_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) 
{ 
	instance* inst = (instance*) userhandle;  
	double* xstar = (double*) malloc(inst->ncols * sizeof(double));
	int *succ = (int*) malloc(inst->nnodes*sizeof(int));
	int ncomp = 0;
	int *comp = (int*) malloc(inst->nnodes*sizeof(int));
	double objval = CPX_INFBOUND; 

	if ( CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols-1, &objval) ) print_error("CPXcallbackgetcandidatepoint error");

	int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread); 
	int mynode = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode); 
	double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent); 
	if ( inst->verbose >= 100 ) printf("\n[Callback] : at node %5d thread %2d incumbent %10.2lf\n", mynode, mythread, incumbent);

	build_sol(xstar, inst, succ, comp, &ncomp);
	
	double curr_cost = compute_sol_cost(inst, succ, inst->nnodes);

	double *xheu = (double *) calloc(inst->ncols, sizeof(double));  // all zeros, initially
	int *ind = (int *) malloc(inst->ncols * sizeof(int));

	if (ncomp > 1){
		for (int k = 1; k <= ncomp; k++){ add_one_sec(inst, NULL, NULL, k, comp, context); }
		if (inst->patching_heur) {
			
			double patch_result = patching_heuristic(inst, NULL, NULL, succ, comp, &ncomp, curr_cost, context);
			if (patch_result < 0) printf("\nThe Patching Heuristic exceeded the time limit, so it did not perform any changes");
			else curr_cost = patch_result;
			if (inst->verbose >= 80) printf("\n[Callback] : cost after patching heuristic is %f\n", curr_cost);
			update_best(inst, succ, curr_cost); // needed in case of a early stop due to low timelimit --> at least we have a tsp feasible solution other than the one given by the heuristic method
			if (curr_cost > 0 && curr_cost < incumbent){
				for ( int j = 0; j < inst->ncols; j++ ) ind[j] = j;
				for ( int i = 0; i < inst->nnodes; i++ ) xheu[xpos(i,succ[i],inst)] = 1.0;
				if ( CPXcallbackpostheursoln(context, inst->ncols, ind, xheu, curr_cost, CPXCALLBACKSOLUTION_NOCHECK) ) print_error("CPXcallbackpostheursoln() error");
				if (inst->verbose >= 80) printf("\n[Callback] : heuristic solution posted\n");
			}
		}
	}
	
	free(ind);
	free(xheu);
	free(comp); 
	free(succ);
	free(xstar);
	return 0; 
}

double Branch_and_Cut(instance *inst, int *best_sol, double timelimit){
	double greedy_timelimit = timelimit/10;
	double BC_timelimit = timelimit/10*9;

	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);

	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;

    double ub = TSPgreedy(inst, best_sol, 0.1, greedy_timelimit, 1);
	double *xheu = (double*) calloc(inst->ncols, sizeof(double));  // all zeros, initially
	for ( int i = 0; i < inst->nnodes; i++ ) {
		xheu[xpos(i, best_sol[i], inst)] = 1.0;
	}
	if (inst->verbose >= 30) printf("\nStarting Upper Bound: %f\n", ub);

	int *ind = (int *) malloc(inst->ncols * sizeof(int));
	for ( int j = 0; j < inst->ncols; j++ ) ind[j] = j;
	int effortlevel = CPX_MIPSTART_NOCHECK;  
	int beg = 0; 			

	if (CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");	
	
	// Cplex's parameter setting
	if (inst->verbose >= 50) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	else CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);

	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);	
	CPXsetdblparam(env, CPX_PARAM_TILIM, BC_timelimit);

	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE; // ... means lazyconstraints
	if ( CPXcallbacksetfunc(env, lp, contextid, incumbent_callback, inst) ) print_error("CPXcallbacksetfunc() error");

	error = CPXmipopt(env,lp);
		if ( error ) 
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}

	if ( CPXgetx(env, lp, xheu, 0, ncols-1) ) printf("CPXgetx() found no better solution\n");	

    CPXgetobjval(env, lp, &ub);
	double best_bound = 0;
	CPXgetbestobjval(env, lp, &best_bound);

	if (inst->verbose >= 30) printf("\nUB: %11.3f   LB: %11.3f\n", ub, best_bound);

	int ncomp = 0;
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	build_sol(xheu, inst, best_sol, comp, &ncomp);

	update_best(inst, best_sol, ub);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	free(comp);
	free(ind);
	free(xheu);

	return ub;
}


// heur_par > 1 means local branching, 0 < heur_par < 1 means hard fixing
double branch_and_cut_HF(instance *inst, int *best_sol, double timelimit, double heur_par){
	double t_start = inst->t_start;
	double t_end = 0;

	double greedy_timelimit = timelimit/10;
	double BC_timelimit = timelimit/10*9;

	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);

	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;

    double ub = TSPgreedy(inst, best_sol, 0.1, greedy_timelimit, 1); 
	double *xheu = (double*) calloc(inst->ncols, sizeof(double));  // all zeros, initially
	for ( int i = 0; i < inst->nnodes; i++ ) {
		xheu[xpos(i, best_sol[i], inst)] = 1.0;
	}
	if (inst->verbose >= 30) printf("\nStarting Upper Bound: %f\n", ub);

	// parameters for adding mip_start
	int *ind = (int *) malloc(inst->ncols * sizeof(int));
	for ( int j = 0; j < inst->ncols; j++ ) ind[j] = j;
	int effortlevel = CPX_MIPSTART_NOCHECK;  
	int beg = 0; 

	int K = 10;

	do {
		if (heur_par > 1){ // local_branching
			add_invalid_cons(inst, env, lp, xheu, K);
		} else if (heur_par > 0){ // hard_fixing
			double one = 1.0;
			char bound = 'L';
			for (int i = 0; i < inst->nnodes; i++){
				if (random01() < heur_par){
					int pos = xpos(i, best_sol[i], inst);
					CPXchgbds(env, lp, 1, &pos, &bound, &one);
				}
			}
		} else {
			print_error("Wrong paramterer value of heur_par");
		}

		if (CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");

		if (inst->verbose >= 50) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
		else CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
		CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);
		double cpx_tl = (BC_timelimit - (second() - t_start) < BC_timelimit/5) ? BC_timelimit - (second() - t_start) : BC_timelimit/5;
		CPXsetdblparam(env, CPX_PARAM_TILIM, cpx_tl);

		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
		if ( CPXcallbacksetfunc(env, lp, contextid, incumbent_callback, inst) ) print_error("CPXcallbacksetfunc() error");

		error = CPXmipopt(env,lp);
		if ( error ) 
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}

		if ( CPXgetx(env, lp, xheu, 0, ncols-1) ) printf("CPXgetx() found no better solution\n");	

		double new_ub;
		CPXgetobjval(env, lp, &new_ub);

		if (heur_par > 1){
			int nrows = CPXgetnumrows(env, lp);
			if (CPXdelrows(env, lp, nrows - 1, nrows - 1)) print_error("Error in deleting local branch");
			if (new_ub >= ub * 0.99999) { // no improvement
				K += 10;
			} else {
				ub = new_ub;
			}
		} else if (heur_par > 0){
			printf("\nRemoving constraints\n");
			for (int i = 0; i < inst->ncols; i++){
				int ind = i;
				double bd = 0.0;
				char lb = 'L';
				if (CPXchgbds(env, lp, 1, &ind, &lb, &bd)) print_error("Error in CPXchgbds()");
			}
			if (new_ub < ub){
				ub = new_ub;
			}	
		}

		int ncomp = 0;
		build_sol(xheu, inst, best_sol, comp, &ncomp);

		update_best(inst, best_sol, compute_sol_cost(inst, best_sol, inst->nnodes));

		t_end = second();

	} while ((t_end - t_start) < BC_timelimit && heur_par > 0);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	free(ind);
	free(xheu);
	free(comp);

	return ub;
}