/*
  This headerfile has my implementation of the LAPACK routines dgesv_ and dgbsv_  for C++.
  
  dgesv:

  This function solves the general set of real linear equations of the
  form Ax = b, where A is a given n by n real matrix and b is a given
  n element vector.  The n element vector x is returned.  

  My function call is of the form

    void dgesv(double **A, double *b, int n, double *x)

    A: the left hand size n by n matrix
    b: the right hand side n element vector
    n: the dimension of A and b
    x: the n element vector for returned values
 
  dgbsv:

  This function solves the system of equation of the form Ax=b where A is
  a sparse banded matrix with dimension of n and b is a given n element vector.
  The n element vector x is returned.

   My function call is of the form

    void dgesv(double **M, double *b, int n, double *x)

    A: nonzero diagonals from updiaonal to low diagonals respectively
    b: the right hand side n element vector
    n: the dimension of A and b
    x: the n element vector for returned values

*/

#ifndef MYLAPACK_HPP
#define MYLAPACK_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

////Lapack Decleration
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda,
                       int *ipiv, double *b, int *ldb, int *info);

extern "C" void dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab,int *ipiv, double *b, int *ldb, int *info);

extern "C" void dgtsv_(int *n, int *nrhs, double *dl, double *d, double*du, double *b, int *ldb, int *info);


class Lapack{

public:

	Lapack(){};

	~Lapack();
     
    void setup(unsigned Num_Channel, unsigned Num_Reservoir, unsigned Total_Cells, int *Channel_Num_Cells, double *dx, int **Connectivity);
    
    void dgesv(double **LHS_2, double *rhs_2, int n);
    
    void dgbsv(double **LHS_1, double *rhs_1, int n, int nrhs_);

    void dgtsv(double *DL, double *D, double *DU, double *B, int n_, int nrhs_);
    
    void  find_three_P_Phi(double *A1_x_, double *B1_x_, double *A2_x_, double *B2_x_, double *f_1_x_, double *f_2_x_, double *S_x_, double *RHS_, double *dx, double *Q1, double *Q2, double *Q3, double *I1, double *I2, double *I3 );
    
    void  find_reservoir_P_Phi(double *Q1, double *Q2, double *Q3, double *I1, double *I2, double *I3, int **Connectivity, int **Connecting_Channel, bool *Reservoir_Pressure_Type, bool *Reservoir_Potential_Type, double *Reservoir_Pressure_Value, double *Reservoir_Potential_Value);
    
    void  find_final_P_Phi_Q_I(double *Q1, double *Q2, double *Q3, double *I1, double *I2, double *I3, int **Channel_End_Reservoir, double *Q, double *I, double *Pressure, double *Potential);

    void P_Phi(double *A1_x_, double *B1_x_, double *A2_x_, double *B2_x_, double *RHS_, double *dx,  int **Channel_End_Reservoir, double *f_1_x_, double *f_2_x_, double *Q, double *I, double *Pressure, double *Potential);
    
    void  solve(double *A1_x_, double *B1_x_, double *A2_x_, double *B2_x_, double *f_1_x_, double *f_2_x_, double *S_x_, double *RHS_, double *dx, int **Connectivity, int **Connecting_Channel, int **Channel_End_Reservoir, bool *Reservoir_Pressure_Type, bool *Reservoir_Potential_Type, double *Reservoir_Pressure_Value, double *Reservoir_Potential_Value, double *Q1, double *Q2, double *Q3, double *I1, double *I2, double *I3, double *Q, double *I, double *Pressure, double *Potential);

    void  find_reservoir_dc(double *Channel_Inlet_Flux_, double *Channel_Outlet_Flux_, double **three_inlet_flux_, double **three_outlet_flux_, bool*Reservoir_Pressure_Type, bool *Reservoir_Potential_Type, int **Connectivity, int **Connecting_Channel, double *Reservoir_Volume, double dt, double *reservoir_dc);
        
private:
   
    int	      *res_connections;
    int       *col; 
    int       Num_Channel_;         //Num of channel
    int       Num_Reservoir_;       //num of reservoirs
    int       *Channel_Num_Cells_;   //num of cells for each channel
    double    *dx_;
    int       Total_Cells_;         //total cells including reservoirs
    int       cell_counter;
    int       face_counter;
    int       kl;                   //num of lower diagonals
    int       ku;                   //num of upper diagonals (for my system kl=ku)
    int       ldab;                 //leading dimension of matrix ab. we take:ldab=ku+2*kl+1
    int       lda;                  //leading dimension of matrix a.
    int       ldb;
    int       nrhs;
    int       size;
    int       info;
    
    double    **LHS_1;              //dimension:ldab*(2*Total_Cells)
    double    *rhs_1;              //dimension: (2*Total_Cells)*3 containts rhs vectors for 3 time solving
    double    *ab;                  //Fortran version of LHS_1 input for dgbsv_
    
    double    **LHS_2;              //dimension: 2*Num_Reservoir* 2*Num_Reservoir
    double    *a;                   //Fortran version of LHS_2 input for dgesv_
    double    *rhs_2;               //rhs for reservoir system of eqn.

    double    **LHS_3;              // dimension: Num_Reservoir* Num_Reservoir 

    int    *ipiv;                //permutation matrix used by lapack
    
    double    **P;                  //dimension: Total_Cells*3 contains P for 3 time solving
    double    **Phi;                //dimension: Total_Cells*3 contains Phi for 3 time solving
    double    *Reservoir_P;         //contains P values in reservoir
    double    *Reservoir_Phi;       //contains Phi values in reservoir

    double    *rhs_p_phi;
    
};
#endif
