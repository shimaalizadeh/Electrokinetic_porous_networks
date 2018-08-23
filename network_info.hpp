#ifndef NETWORK_INFO_HPP
#define NETWORK_INFO_HPP

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<string>
#include<stdio.h>
#include<assert.h>
using namespace std;

class Network{

 public:

   Network(){};
  ~Network();
  
   void setup(int Num_Channel, int Num_Reservoir, int restart, double dt, double Tmax, double time, double period, string network_input_file);

    void Initial_Memory_Allocation(int Num_Channel,int Num_Reservoir);
    void Memory_Allocation();
    void Read_filename_1(string filename_1);
    void read_restart_files(void);
    void set_electrochemical_parameters();
    void set_to_zero(double *vector, int size);

    void set_dt(void);
    void compute_gradients(void);
    void store_reservoir_data(double *data, double*Reservoir_data);
    void update_C_bar(int myid, double tt);
    void update_Old_Cbar(void);
    void compute_first_iteration_rhs_cbar_part(void);
    void compute_rhs_cbar_part(void);
    void write_restart_network(double time);
    void write_restart_channel(double time);
    double check_C_conservativity(void);

   // for initial file reading
   double*  sigma_val;
   double*  lambda_val;
   double*  s_val;
   double*  c0_val;
  
   //public variables
   double       slit_gp_bar;          //Average Pressure driven flow in a slit
   double       circle_gp_bar;
   double       Pe;
   double       lambda_ref;
   int          Num_Channel;             
   int          Num_Reservoir;
   int          Total_Cells;             //including reservoir cells
   int          Total_Faces;
   int         *Channel_type;
   int         *Channel_blocked;
   int         *Reservoir_type;
   int         *Channel_Num_Cells;       //excluding reservoir cells 
   int         **Channel_End_Reservoir;  //num_channel*2 matrix: 1st column is inlet reservoir, 2nd column is end reservoir
  int         **Connectivity;
  int         **Connecting_Channel;
  bool        *Reservoir_Pressure_Type; //if 0--> end reservoir if 1-->internal
  bool        *Reservoir_Potential_Type;//if 0--> end reservoir if 1-->internal
  double      *Reservoir_Pressure_Value;
  double      *Reservoir_Potential_Value;
  double      *Reservoir_C_bar;
  double      *Reservoir_C0;
  double      *Reservoir_Sigma;
  double      *Reservoir_Lambda;
  double      *Reservoir_Pressure;
  double      *Reservoir_Potential;
  double      *Reservoir_Volume;
  double      *dx;
  double      *C_bar;
  double      *C0;
  double      *Cs;
  double      *C_bar_n;  //cbar in moment (n)
  double      *C_bar_n_1; //cbar in moment (n-1)
  double      *rhs_cbar_part;
  double      *sigma_star;
  double      *lambda;
  double      *lambda_for_table;
  double      *nondimensional_Sp;
  double      *Pressure;
  double      *Potential;
  double      *U;
  double      *Q1, *Q2, *Q3;
  double      *I1, *I2, *I3;
  double      *Q, *I;

  double      Tmax;
  double      time;
  double      dt;
  int         period;

  double      *dC_bar;

  //for simplicity in flux calculation
  double      *dPdx;
  double      *dMudx;
 
 private: 
  int         restart;
  int         cell_counter;
  int         face_counter;
  
  //electrochemical parameters
  int    z;                 //valence
  double e;                 //elementary charge(Coulombs!)
  double Kb;                //Boltzmann Constant
  double epsilon0;          //permitivity constant
  double epsilon;           //dielectric relative permitivity
  double mu;                //dynamic viscosity(Pa.s)
  double T;                 //temperature(K)
  double Dif;               //mass Diffusivity (m^2/s),
  double C_ref;             //reference concentration(Molar)
  double NA;                //Avogadro number
  
  };
#endif 
