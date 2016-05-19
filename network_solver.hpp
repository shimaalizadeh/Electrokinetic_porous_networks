#ifndef NETWORK_SOLVER_HPP
#define NETWORK_SOLVER_HPP

#include "network_info.hpp"
#include "table.hpp"
#include "mylapack.hpp"

class NetworkSolver{

	public:
	
	NetworkSolver(int myid, int Num_Channel, int Num_Reservoir, int restart, double dt, double Tmax, double time, int period, double right_voltage);
	
	~NetworkSolver();
	
	void reset(void);
	
	void solve(void);

        void report_status();
	
	private:
    int myid;
	int Num_Channel;
	int Num_Reservoir;
	int restart;
	double dt;
	double Tmax;
	double time;
	int period;
	double right_voltage;
	
	Network network;
	Table table;
	Lapack solver;

};
#endif
