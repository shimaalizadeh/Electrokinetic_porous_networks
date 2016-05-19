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
NetworkSolver::NetworkSolver(int myid, int Num_Channel, int Num_Reservoir, int restart, double dt, 
double Tmax, double time, int period, double right_voltage)
{
	this-> myid = myid;
	this-> Num_Channel = Num_Channel;
	this-> Num_Reservoir = Num_Reservoir;
	this-> restart = restart;
	this-> period = period;
	this-> time = time;
	this-> dt = dt;
	this-> Tmax = Tmax;
	this-> right_voltage = right_voltage;
	
	////1. network object
  	network.setup("input_files/Network.in", Num_Channel, Num_Reservoir, restart, right_voltage);
  	network.set_simulation_time(time, dt, Tmax, period);

  	////2.table reading object
  	table.setup(Num_Channel, Num_Reservoir, network.Channel_Num_Cells, network.Channel_type,
  	 network.Reservoir_type, network.lambda_ref, network.Pe, network.slit_gp_bar, 
  	 network.circle_gp_bar, "/home/shima86/tables/slit_sigma_lambda0.in", "/home/shima86/tables/slit_sigma_lambdastar.in",
  	 "/home/shima86/tables/circle_sigma_lambda0.in", "/home/shima86/tables/circle_sigma_lambdastar.in");

  	////3. poisson solver object
  	solver.setup(Num_Channel, Num_Reservoir, network.Total_Cells, network.Channel_Num_Cells, 
		network.dx);

	Mesh_Generation(Num_Channel, Num_Reservoir, network.Ns, network.Nt, network.Nz, network.ds, network.dn, network.dz, network.theta, network.x0, network.y0, network.z0);	
}

//destructor
NetworkSolver::~NetworkSolver(){}


void NetworkSolver::reset(void)
{
	restart = 0;
	time = 0.;
	
	////1. network object
  	network.reset("input_files/Network.in", right_voltage);
  	network.set_simulation_time(time, dt, Tmax, period);

  	return;
}

void NetworkSolver::solve()
{
  
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
  	network.store_reservoir_data(network.C0, network.Reservoir_C0);
  	network.store_reservoir_data(network.Pressure, network.Reservoir_Pressure);
  	network.store_reservoir_data(network.Potential, network.Reservoir_Potential);


  	////Initial writing
  	int counter = 0;
  	write(Num_Channel, Num_Reservoir, network.Ns, network.Nt, network.Nz, 
			network.C_bar, network.Reservoir_C_bar, network.C0, network.Reservoir_C0, 
			network.Pressure, network.Reservoir_Pressure, network.Potential, 
			network.Reservoir_Potential, counter);
  				    
  while(time< Tmax){
	
		time += dt;

       // shock treatment based on C_new and I_old
       network.Shock_Treatment();

       table.do_table_reading(network.sigma_star,network.lambda_for_table, network.lambda, 
			      network.C_bar, network.C0, network.Reservoir_C0, 
			      network.Channel_End_Reservoir, network.Cs, network.nondimensional_Sp,
			      network.dx, network.Art_Diff, time);
			      
       solver.solve(table.A1_x_, table.B1_x_, table.A2_x_, table.B2_x_, table.f1_x_, table.f2_x_, 
		    table.S_x_, table.RHS_, network.dx, network.Connectivity, 
		    network.Connecting_Channel, network.Channel_End_Reservoir, 
		    network.Reservoir_Pressure_Type, network.Reservoir_Potential_Type, 
		    network.Reservoir_Pressure_Value, network.Reservoir_Potential_Value, 
		    network.Q1, network.Q2, network.Q3, network.I1, network.I2, network.I3, 
		    network.Q, network.I, network.Pressure, network.Potential);

		if(int(time / dt) % period == 0){
			
			counter++;
			
			write(Num_Channel, Num_Reservoir, network.Ns, network.Nt, network.Nz, 
			network.C_bar, network.Reservoir_C_bar, network.C0, network.Reservoir_C0, 
			network.Pressure, network.Reservoir_Pressure, network.Potential, 
			network.Reservoir_Potential, counter);
		}
       
       // update reservoir pressure and potential
       network.store_reservoir_data(network.Pressure, network.Reservoir_Pressure);
       network.store_reservoir_data(network.Potential, network.Reservoir_Potential);
	
       // compute dPdx and dMudx
       network.compute_gradients();
       
       // shock treatment based on C_new and I_new
       network.Shock_Treatment();

       table.find_flux_n(network.C_bar, network.C0, network.dPdx, network.dMudx,
			 network.nondimensional_Sp, network.dx, network.U, time);

       table.solve_for_dc(solver, network.C_bar, network.dMudx, network.U, network.nondimensional_Sp, 
			  network.Reservoir_Pressure_Type, network.Reservoir_Potential_Type, 
			  network.Reservoir_Volume, network.Connectivity, network.Connecting_Channel, 
			  network.Channel_End_Reservoir, network.dx, dt, network.dC_bar);
   
       network.Update_C_bar(myid, time);

       report_status();
 
  }//End of while

	////Final writing
    counter++;
	write(Num_Channel, Num_Reservoir, network.Ns, network.Nt, network.Nz, 
			network.C_bar, network.Reservoir_C_bar, network.C0, network.Reservoir_C0, 
			network.Pressure, network.Reservoir_Pressure, network.Potential, 
			network.Reservoir_Potential, counter);
		
    return;
 
}

void  NetworkSolver::report_status(){
	
	cout<<"time= "<<time<<endl;
	cout<<"I[0] = "<<network.I[2]<<"  , I[1]= "<<network.I[3]<<"    , I[3]= "<<network.I[5]<<"I[end]= "<<network.I[16]<<endl;
	cout<<"Q[0] = "<<network.Q[2]<<"  , Q[1]= "<<network.Q[3]<<"  , Q[2]= "<<network.Q[5]<<"  , Q[3]= "<<network.Q[9]<<"    , Q[end]= "<<network.Q[16]<<endl;
	cout<<endl;

	return;
}
