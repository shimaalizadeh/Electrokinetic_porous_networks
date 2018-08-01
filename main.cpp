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
    if(argc != 9){
	cerr<<"ERROR: 13 arguments required: \n"<<
    "# of pores, \n"<<
    "# of reservoirs, \n"<<
    "# of horizontal pores in x, \n"<<
    "# of horizontal pores in y, \n"<<
    "total # of horizontal pores, \n"<<
    "delta_t, \n"<<
    "max_time, \n"<<
    "write_frequency, \n"<<
    "restart(0 or 1)"<<
    "table_path\n"<<endl;
        exit(EXIT_FAILURE);
    }

    // some default values
    double current_time = 0.;
    double right_voltage = 0.; 

    // get arguments
    int num_pores = atoi(argv[1]);
    int num_reservoirs = atoi(argv[2]);
    int hp_x = atoi(argv[3]);
    int hp_y = atoi(argv[4]);
    int h_pores = atoi(argv[5]);
    double dt = atof(argv[6]);
    double tmax = atof(argv[7]);
    int period = atoi(argv[8]);
    int restart = atoi(argv[9]);
    string table_path = argv[10];
    
    //setup the solver
    NetworkSolver network_solver(myid, 
                                 num_pores, 
                                 num_reservoirs, 
                                 hp_x, hp_y, h_pores, 
                                 restart, dt, tmax, 
                                 current_time, period, 
                                 right_voltage,
                                 table_path);
 
   // solve for the first time step
   network_solver.solve_first_iteration();
     
   // solve for other timesteps
   network_solver.solve();

   return(0);
}
