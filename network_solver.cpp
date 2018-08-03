#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<stdio.h>
#include "network_solver.hpp"
#include "network_write.hpp"

using namespace std;

//Constructor
NetworkSolver::NetworkSolver(int myid, int Num_Channel, int Num_Reservoir, int hp_x, int hp_y, int h_pores, int restart, double dt, double Tmax, double time, int period, double right_voltage, string table_path){

	this-> myid = myid;
	this-> Num_Channel = Num_Channel;
	this-> Num_Reservoir = Num_Reservoir;
	this-> hp_x = hp_x;
	this-> hp_y = hp_y;
	this-> h_pores = h_pores;
	this-> restart = restart;
	this-> period = period;
	this-> dt = dt;
	this-> Tmax = Tmax;
	this->time = time;
	this-> right_voltage = right_voltage;

	restart_period = 4*period;
	inner_iteration = 3;

	////1. network object
  	network.setup(Num_Channel, Num_Reservoir, restart, dt, Tmax, time, period, right_voltage);
	this->time = network.time;

	if( int(this->time/dt) % period == 0){
        	counter = int(this->time/(period*dt)); // used for writing
	}
	else{
		counter = int(this->time/(period*dt))+1;
	}

  	////2.table reading object
  	table.setup(Num_Channel, Num_Reservoir, h_pores, hp_x, hp_y, network.Channel_Num_Cells, network.Channel_blocked, network.Channel_type,
  	 network.Reservoir_type, network.lambda_ref, network.Pe, network.slit_gp_bar, 
  	 network.circle_gp_bar, table_path+"slit_sigma_lambda0.in", table_path+"slit_sigma_lambdastar.in",

  	 table_path+"circle_sigma_lambda0.in", table_path+"circle_sigma_lambdastar.in");

  	////3. poisson solver object
  	solver.setup(Num_Channel, Num_Reservoir, network.Total_Cells, network.Channel_Num_Cells, 
		network.dx, network.Connectivity);

}

//destructor
NetworkSolver::~NetworkSolver(){}

void NetworkSolver::solve()
{
    time_scheme_factor = 3. / 2.; //for initial time step it's 1, later it's 3./2.

  	while(time< Tmax){

		time += dt;

        	//// start inner loop
		network.update_Old_Cbar();

		for(int i=0; i<inner_iteration; i++){

			network.compute_rhs_cbar_part();

			table.do_table_reading(network.sigma_star,network.lambda_for_table, network.lambda, 

			network.C_bar, network.C0, network.Reservoir_C0, network.Channel_End_Reservoir, 

			network.Cs, network.nondimensional_Sp, network.dx, time);

			solver.solve(table.A1_x_, table.B1_x_, table.A2_x_, table.B2_x_, table.f1_x_, 

			table.f2_x_, table.S_x_, table.RHS_, network.dx, network.Connectivity, 

			network.Connecting_Channel, network.Channel_End_Reservoir, 

			network.Reservoir_Pressure_Type, network.Reservoir_Potential_Type, 

			network.Reservoir_Pressure_Value, network.Reservoir_Potential_Value, 

			network.Q1, network.Q2, network.Q3, network.I1, network.I2, network.I3, 

			network.Q, network.I, network.Pressure, network.Potential);

       			// update reservoir pressure and potential
       			network.store_reservoir_data(network.Pressure, network.Reservoir_Pressure);
       			network.store_reservoir_data(network.Potential, network.Reservoir_Potential);
	
       			// compute dPdx and dMudx
       			network.compute_gradients();

       			table.find_flux_n(network.C_bar, network.C0, network.dPdx, network.dMudx, 
			network.nondimensional_Sp, network.dx, network.U, time);

			table.solve_for_dc(solver, network.C_bar, network.rhs_cbar_part, network.dMudx, 
			network.U, network.nondimensional_Sp, network.Reservoir_Pressure_Type, 
			network.Reservoir_Potential_Type, network.Reservoir_Volume, network.Connectivity, 
			network.Connecting_Channel, network.Channel_End_Reservoir, network.dx, dt,
			time_scheme_factor, network.dC_bar);

       			network.update_C_bar(myid, time);

			////Storing Reservoir C_bar
			network.store_reservoir_data(network.C_bar, network.Reservoir_C_bar);

		}
#if 0
		if( int(time/dt) % network.period == 0){

			write(Num_Channel, Num_Reservoir, network.Ns, network.Nt, network.Nz, 
			network.C_bar, network.Reservoir_C_bar, network.C0, network.Reservoir_C0, 
			network.Pressure, network.Reservoir_Pressure, network.Potential,
			network.Reservoir_Potential, network.U, h_pores, counter);

			counter++;

   	        }
#endif
		if( int(time/dt) % restart_period == 0){

			//write restart files
			cout<<"write restart files at time "<<time<<endl;
                        network.write_restart_network(time);
                        network.write_restart_channel(time);
                }


	

       		report_status();

  	}//End of while

  	////Final writing
	write(Num_Channel, Num_Reservoir, network.Ns, network.Nt, network.Nz, network.C_bar, 
	network.Reservoir_C_bar, network.C0, network.Reservoir_C0, network.Pressure, 
	network.Reservoir_Pressure, network.Potential, network.Reservoir_Potential, network.U, h_pores, counter);

	////write currents, flowrates, area cross-sections
	write_current_flow(Num_Channel, network.Q, network.I, network.nondimensional_Sp);

	////Write restart files
	network.write_restart_network(time);
	network.write_restart_channel(time);
	cout<<"restart files written!"<<endl;

	return; 
}

void NetworkSolver::solve_first_iteration(void){

	cout<<"running the first iteration..."<<endl;
  	time_scheme_factor = 1.;


	////find initial C_bar if restart==0
  	if(restart == 0)
  	{

		table.Find_Initial_C_bar(network.sigma_star, network.lambda_for_table, network.C0,
			         network.Channel_End_Reservoir, network.C_bar, 
				 network.Reservoir_C_bar, network.time);

  	}
  	else
  	{
		table.Import_Table_2(network.lambda_for_table);
  	}

  	////store reservoir data
  	network.store_reservoir_data(network.C_bar, network.Reservoir_C_bar);
  	network.store_reservoir_data(network.Pressure, network.Reservoir_Pressure);
  	network.store_reservoir_data(network.Potential, network.Reservoir_Potential);
#if 0
	////initial writing
	if( int(time/dt) % network.period == 0){
		write(Num_Channel, Num_Reservoir, network.Ns, network.Nt, network.Nz, network.C_bar, 
		network.Reservoir_C_bar, network.C0, network.Reservoir_C0, network.Pressure, 
		network.Reservoir_Pressure, network.Potential, network.Reservoir_Potential, network.U, h_pores, counter);
		counter++;
	}
#endif
        //// start inner loop
	time += dt;
	network.update_Old_Cbar();
	for(int i=0; i<inner_iteration; i++){

		network.compute_first_iteration_rhs_cbar_part();

		table.do_table_reading(network.sigma_star,network.lambda_for_table, network.lambda,
		network.C_bar, network.C0, network.Reservoir_C0, network.Channel_End_Reservoir, 
		network.Cs, network.nondimensional_Sp, network.dx, time);

		solver.solve(table.A1_x_, table.B1_x_, table.A2_x_, table.B2_x_, table.f1_x_, 
		table.f2_x_, table.S_x_, table.RHS_, network.dx, network.Connectivity, 
		network.Connecting_Channel, network.Channel_End_Reservoir, network.Reservoir_Pressure_Type,
		network.Reservoir_Potential_Type, network.Reservoir_Pressure_Value, 
		network.Reservoir_Potential_Value, network.Q1, network.Q2, network.Q3, network.I1,
		network.I2, network.I3, network.Q, network.I, network.Pressure, network.Potential);

       		// update reservoir pressure and potential
       		network.store_reservoir_data(network.Pressure, network.Reservoir_Pressure);
       		network.store_reservoir_data(network.Potential, network.Reservoir_Potential);

       		// compute dPdx and dMudx
       		network.compute_gradients();

       		table.find_flux_n(network.C_bar, network.C0, network.dPdx, network.dMudx, 
		network.nondimensional_Sp, network.dx, network.U, time);

		table.solve_for_dc(solver, network.C_bar, network.rhs_cbar_part, network.dMudx, 
		network.U, network.nondimensional_Sp, network.Reservoir_Pressure_Type, 
		network.Reservoir_Potential_Type, network.Reservoir_Volume, network.Connectivity, 
		network.Connecting_Channel, network.Channel_End_Reservoir, network.dx, dt, 
		time_scheme_factor, network.dC_bar);

		network.update_C_bar(myid, time);
				
		////Storing Reservoir C_bar
		network.store_reservoir_data(network.C_bar, network.Reservoir_C_bar);
	}

       	report_status();
  return;

}

void  NetworkSolver::report_status(){

	cout<<"time= "<<time<<endl;
    double tmp = 0.0;
	for(int i = 0; i < Num_Channel; i++)
		tmp += network.I[i];

		cout<<"I_tot = "<<tmp<<endl;
		cout<<endl;
	return;
}
