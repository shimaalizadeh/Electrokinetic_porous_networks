#ifndef NETWORK_INFO_HPP
#define NETWORK_INFO_HPP

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<string>
#include<stdio.h>
using namespace std;

#include"mymath_tools.h"

class Network{

 public:

   Network(){};
  ~Network();
  
   void setup(const char* filename_1, int Num_Channel_, int Num_reservoir_, int restart_, double right_voltage);

    void reset(const char* filename_1, double right_voltage);
    
    void reset_from_file(const char* filename_1, double right_voltage);
    
    void Initial_Memory_Allocation(int Num_Channel,int Num_Reservoir);
    void Read_filename_1(const char* filename_1, double right_voltage);
    void Memory_Allocation();
    void set_electrochemical_parameters();
    void set_to_zero(double *vector, int size);

    void set_dt();
    void set_simulation_time(double time, double dt, double T_max, int period);
    void store_reservoir_data(double *data, double*Reservoir_data);
    void Update_C_bar(int myid, double tt);
    double* Shock_Treatment();      
    void Write_filename_1(const char* filename_1);
    double check_C_conservativity(void);
    void   compute_gradients(void);

   // for initial file reading
   double* sigma_val;
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
  double      *sigma_star;
  double      *lambda;
  double      *lambda_for_table;
  double      *nondimensional_Sp;
  double      *Pressure, *Potential;
  double      *U;
  double      *Q1, *Q2, *Q3;
  double      *I1, *I2, *I3;
  double      *Q, *I;

  double      T_max;
  double      time;
  double      dt;
  int         period;

  double      *dC_bar;
  double      *Art_Diff;     //artificial diffusivity

  //for simplicity in flux calculation
  double      *dPdx;
  double      *dMudx;
 
  //geometry
  int         *Ns;                      //# of cells in streamwise
  int         *Nt;                      //# of cells in transverse
  int         *Nz;                      //# of cells in z
  double      *ds;
  double      *dn;
  double      *dz;
  double      *x0;
  double      *y0;
  double      *z0;
  double      *theta;

 private: 
  int         restart;
  int         cell_counter;
  int         face_counter;
  int         Art_Diff_Switch; 
  
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
