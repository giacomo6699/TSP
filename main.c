#include "tsp_solver.h"
#include <cplex.h>

double second();
void print_error(const char *err);
int main(int argc, char **argv);
void read_input(instance *inst);
void get_random_inst(instance *inst);    
double random01();     
void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst); 
void print_instance(instance *inst);
void print_travel(int *succ, instance *inst);

double get_random_solution(instance *inst, int *succ);
void print_succ(int *succ, int n);

void debug(const char *err) { printf("\nDEBUG: %s \n", err); fflush(NULL); }

int number_of_nonempty_lines(const char *file)  // warning: the last line NOT counted if it is does not terminate with \n (as it happens with some editors) 
{
	FILE *fin = fopen(file, "r");
	if ( fin == NULL ) return 0;
	char line[123456]; 
	int count = 0;
	while( fgets(line, sizeof(line), fin) != NULL ) { printf(" len %4d\n", (int) strlen(line)); if ( strlen(line) > 1 ) count++; }
	fclose(fin);   
	return count; 	
}


void free_instance(instance *inst)
{     
	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->best_sol);
}

void get_random_inst(instance *inst){
	if (inst->randomseed != 0){
		srand(82356+abs(inst->randomseed));
		for (int k = 0; k < 1000; k++) rand();
	}
	inst->xcoord = malloc(inst->nnodes * sizeof(double));
	inst->ycoord = malloc(inst->nnodes * sizeof(double));
	for (int i = 0; i < inst->nnodes; i++){
		inst->xcoord[i] = ((double) rand() / RAND_MAX) * 1000;
		inst->ycoord[i] = ((double) rand() / RAND_MAX) * 1000;
	}
}    

int main(int argc, char **argv) 
{
	instance inst;
	inst.t_start = second();

	parse_command_line(argc,argv, &inst);   

	if (inst.nnodes >= 0) 
		get_random_inst(&inst);
	else 
		read_input(&inst);
	
	int *best_sol = (int*) malloc(inst.nnodes * sizeof(int));

	double benders = benders_loop(&inst, best_sol, 30);
	printf("TSP_n%d_seed%d , %f", inst.nnodes, inst.randomseed, benders);


    
	if (VERBOSE >= 10)   
	{
		printf("... TSP problem solved in %lf sec.s with solution\n", second()-inst.t_start, benders);  
	}
	
	free(best_sol);
	free_instance(&inst);
	return 0; 
}         

void read_input(instance *inst) // simplified CVRP parser, not all SECTIONs detected  
{
                            
	FILE *fin = fopen(inst->input_file, "r");
	if ( fin == NULL ) print_error(" input file not found!");
	
	inst->nnodes = -1;

	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 
	
	int do_print = ( VERBOSE >= 1000 );

	while ( fgets(line, sizeof(line), fin) != NULL ) 
	{
		if ( VERBOSE >= 2000 ) { printf("%s",line); fflush(NULL); }
		if ( strlen(line) <= 1 ) continue; // skip empty lines
	    par_name = strtok(line, " :");
		if ( VERBOSE >= 3000 ) { printf("parameter \"%s\" ",par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 ) 
		{
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "COMMENT", 7) == 0 ) 
		{
			active_section = 0;   
			token1 = strtok(NULL, "");  
			continue;
		}   
		
		if ( strncmp(par_name, "TYPE", 4) == 0 ) 
		{
			token1 = strtok(NULL, " :");  
			if ( strncmp(token1, "CVRP",4) != 0 ) print_error(" format error:  only TYPE == CVRP implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}
		

		if ( strncmp(par_name, "DIMENSION", 9) == 0 ) 
		{
			if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes);	 
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));    
			active_section = 0;  
			continue;
		}


		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 ) 
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "EUC_2D", 6) != 0 ) print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}            
		
		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;   
			continue;
		}
		
		if ( strncmp(par_name, "EOF", 3) == 0 ) 
		{
			active_section = 0;
			break;
		}
		
			
		if ( active_section == 1 ) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");     
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]); 
			continue;
		}     
		
		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the current simplified parser!!!!!!!!!");     
		    
	}                

	fclose(fin);    
	
}

void parse_command_line(int argc, char** argv, instance *inst) 
{ 
	if (argc < 2){
		printf("Usage: %s -help for help\n", argv[0]);
	}
	if (VERBOSE >= 2){
		for (int a = 0; a < argc; a++){
			printf("%s ", argv[a]);
		}
		printf("\n");
	}
	if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1); 
		
	// default   
	strcpy(inst->input_file, "NULL");
	inst->randomseed = 0; 
	inst->timelimit = INFINITY;
	inst->nnodes = 5;
	inst->best_value = INFINITY;
	inst->grasp_par = 1;
	inst->verbose = 10;
	inst->integer_costs = 0;
	inst->patching_heur = 0;
    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-tl") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }				// total time limit
		if ( strcmp(argv[i],"-n") == 0 ) { inst->nnodes = atoi(argv[++i]); continue; } 					// nnodes
		if ( strcmp(argv[i],"-grasp") == 0 ) { inst->grasp_par = atoi(argv[++i]); continue; } 		
		if ( strcmp(argv[i],"-patch_heur") == 0 ) { inst->patching_heur = atoi(argv[++i]); continue; } 		
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-verbose") == 0 ) { inst->verbose = abs(atoi(argv[++i])); continue; }
		if ( strcmp(argv[i],"-int_costs") == 0 ) { inst->integer_costs = 1; continue; }
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }  

	inst->best_sol = (int*) malloc(inst->nnodes*sizeof(int));    

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("\n\nAvailable parameters (0.0) --------------------------------------------------------------------\n");
		printf("-file %s\n", inst->input_file); 
		printf("-tl %lf\n", inst->timelimit); 
		printf("-seed %d\n", inst->randomseed); 
		printf("-n %d\n", inst->nnodes); 
		printf("-verbose %d\n", inst->verbose); 
		printf("-int_costs %d\n", inst->integer_costs); 
		printf("\nEnter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}        
	
	if ( help ) exit(1);

}    




