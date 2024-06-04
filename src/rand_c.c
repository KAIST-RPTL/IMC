#include <stdio.h>
#include <omp.h>
#include <time.h>

int seed = 77777777;
#pragma omp threadprivate(seed)

double rang() {
	
	time_t current_time;
    current_time = time(NULL);
	
	//seed = seed+125346*(int)current_time;
	//seed = seed+125346*omp_get_thread_num(); 
	seed = seed+125346*omp_get_thread_num() + (int)current_time;
	// Initialise the random number generator with different seed in each thread
	// The following constants are chosen arbitrarily... use something more sensible
	//#pragma omp parallel firstprivate(seed)
	//seed = 25234 + 17*omp_get_thread_num();
	
	double x = (double)rand_r(&seed) / (double) RAND_MAX;
	return x; 
}
