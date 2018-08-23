#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include "network_solver.hpp"
using namespace std;


int main(int argc, char *argv[]){
    
    //////////////////////////Initialize MPI///////////////////////
    int myid=0, num_procs=0;   
  
    // get arguments
    if(argc != 10){
	cerr<<"ERROR: 9 arguments required: \n"<<
    "# of pores, \n"<<
    "# of reservoirs, \n"<<
    "total # of horizontal pores, \n"<<
    "delta_t, \n"<<
    "max_time, \n"<<
    "write_frequency, \n"<<
    "restart(0 or 1)"<<
    "table_path\n"<<
    "network_input_file\n"<<endl;
        exit(EXIT_FAILURE);
    }

    // some default values
    double current_time = 0.;

    // get arguments
    int num_pores = atoi(argv[1]);
    int num_reservoirs = atoi(argv[2]);
    int h_pores = atoi(argv[3]);
    double dt = atof(argv[4]);
    double tmax = atof(argv[5]);
    int period = atoi(argv[6]);
    int restart = atoi(argv[7]);
    string table_path = argv[8];
    string network_input_file = argv[9];
    
    //setup the solver
    NetworkSolver network_solver(myid, 
                                 num_pores, 
                                 num_reservoirs, 
                                 h_pores, 
                                 restart, dt, tmax, 
                                 current_time, period, 
                                 table_path,
                                 network_input_file);
 
   // solve for the first time step
   network_solver.solve_first_iteration();
     
   // solve for other timesteps
   network_solver.solve();

   return(0);
}
