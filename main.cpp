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
#include "config.hpp"
using namespace std;


int main(int argc, char *argv[]){
    
    //////////////////////////Initialize MPI///////////////////////
    int myid=0;
  
    // load run settings: either a TOML file or the legacy 9 positional args
    if(argc != 2 && argc != 10){
	cerr<<"Usage:\n"<<
    "  ek_solver run.toml\n"<<
    "  ek_solver <# pores> <# reservoirs> <# horizontal pores> <delta_t> <max_time>\n"<<
    "            <write_frequency> <restart 0|1> <table_path> <network_input_file>\n"<<endl;
        exit(EXIT_FAILURE);
    }

    Config cfg;
    string network_input_file;
    string table_path;
    try{
        cfg = (argc == 2) ? load_config_toml(argv[1]) : load_config_args(argc, argv);
        validate_config(cfg);
        network_input_file = cfg.network_name_for_setup();
        table_path = cfg.tables_prefix();
    }catch(const std::exception& e){
        cerr<<"ERROR: "<<e.what()<<endl;
        exit(EXIT_FAILURE);
    }

    // some default values
    double current_time = 0.;

    int num_pores = cfg.num_pores;
    int num_reservoirs = cfg.num_reservoirs;
    int h_pores = cfg.num_horizontal_pores;
    double dt = cfg.dt;
    double tmax = cfg.max_time;
    int period = cfg.write_every;
    int restart = cfg.restart ? 1 : 0;
    
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
