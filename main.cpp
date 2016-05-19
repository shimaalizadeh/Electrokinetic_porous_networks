#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>
#include "network_solver.hpp"
using namespace std;


int main(int argc, char *argv[])
{
    
    //////////////////////////Initialize MPI///////////////////////
    
    int myid, num_procs;   
  
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   
    // read input files
    ifstream reader("input_files/RightVoltage.in");
    
    int n_voltage;
    
    reader>>n_voltage;
    assert(n_voltage == num_procs);
    int local_sim = 1;
    
    double* voltage = new double [n_voltage];
    
    for(int i=0 ; i<n_voltage ; i++) {
    
    	reader>>voltage[i];
    	
    }
    
    //set-up the solver
    NetworkSolver network_solver(myid, 17, 12, 0, 5e-5, 2., 0., 100, voltage[myid]);
    
    
    // solve h-junction
    network_solver.solve();
    	
    //free memory
    delete [] voltage;
    
    MPI_Finalize();
    return(0);
}
    		
    
