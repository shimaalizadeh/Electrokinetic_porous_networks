#ifndef NETWORK_SOLVER_HPP
#define NETWORK_SOLVER_HPP

#include<string>
#include "network_info.hpp"
#include "table.hpp"
#include "mylapack.hpp"


class NetworkSolver{

	public:
	
	NetworkSolver(int myid, int Num_Channel, int Num_Reservoir, int h_pores, int restart, double dt, double Tmax, double time, int period, string table_path, string network_input_file);
	
	~NetworkSolver();
	
	void solve(void);

	void solve_first_iteration(void);

        void report_status();
	
	private:
	int myid;
	int Num_Channel;
	int Num_Reservoir;
	int h_pores;
	int restart;
	double dt;
	double Tmax;
	double time;
	int period;
	int restart_period;
	double right_voltage;
	int counter;
	int inner_iteration;
	double time_scheme_factor;
	
	Network network;
	Table table;
	Lapack solver;

};
#endif
