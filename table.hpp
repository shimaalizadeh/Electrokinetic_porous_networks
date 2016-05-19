#ifndef TABLE_HPP
#define TABLE_HPP

#include "mylapack.hpp"
using namespace std;

class Table{

 public:
 
   Table(){};
 
  ~Table();

  void setup(int Num_Channel, int Num_Reservoir, int *Channel_Num_Cells, int *Channel_type, int *Reservoir_type, double lambda_ref=1e-9, double Pe=0.7, double slit_gp=-1./3, double circle_gp=-1./2, const char* filename_1="slit_sigma_lambda0.in", const char* filename_2="slit_sigma_lambda_star.in", const char* filename_3="circle_sigma_lambda0.in", const char* filename_4="circle_sigma_lambda_star.in", unsigned max_sigma_star_counter=48, unsigned max_lambda_star_counter=38, double sigma_star_min=-0.1, double lambda_star_min=0.01);

  void Import_Table_1(double *lambda_for_table);

  void Import_Table_2(double *lambda_for_table);

  void Compute_asymptotic(double *lambda_for_table);

  void Find_From_Table(double sigma_star, double lambda, double C_bar, int type, int idx);
    
  void Find_Initial_C_bar(double *sigma_star, double *lambda_for_table, double *C0, int **Channel_End_Reservoir, double *C_bar, double *Reservoir_C_bar, double time);

  void do_table_reading(double *sigma_star,double *lambda_for_table, double *lambda, double *C_bar, double *C0, double *Reservoir_C0,int **Channel_End_Reservoir, double *Cs, double *nondimensional_Sp, double *dx, double *Art_Diff, double time);
    
  void calculate_A_B_f_RHS(double *C_bar, double *C0, double *Cs, double *lambda, double *nondimensional_Sp, double *dx, double *Art_Diff, int **Channel_End_Reservoir);

  void find_flux_n(double *C_bar, double *C0, double *dPdx, double *dMudx,double *nondimensional_Sp, double *dx, double *U, double time);

  void solve_for_dc(Lapack  &solver, double *C_bar, double *dMudx, double *U, double *nondimensional_Sp, bool *Reservoir_Pressure_Type,  bool *Reservoir_Potential_Type, double *Reservoir_Volume, int **Connectivity, int **Connecting_Channel, int **Channel_End_Reservoir, double *dx, double dt, double *dC_bar);

  const double * give_f_bar(void){return(f_bar_);}

    int *Channel_type_;

    int *Reservoir_type_;

    ////some vectors required for future calculations in faces
    double * A1_x_;
    double * A2_x_;
    double * B1_x_;
    double * B2_x_;
    double * f1_x_;
    double * f2_x_;
    double * RHS_;
    double * S_x_;    
    double * sqrt_Cf_x_;
    double * dCfdx;
    private:
    double * A1_;
    double * A2_;
    double * B1_;
    double * B2_;
    double * f1_1_;
    double * f2_1_;

    double * f_bar_;
    double * ge_bar_;
    double * gc_bar_;
    double * gp_m_;
    double * ge_m_;
    double * gc_m_;
    double * gp_p_;
    double * ge_p_;
    double * gc_p_;

    // extra arrays
    double *gcbar_over_fbar_;
    double *gcp_over_fbar_;
    double *gcm_over_fbar_;
    double *sqrt_C_over_fbar_;
    double *fbar_sqrt_C_;

    double * Gp_m_;
    double * Ge_m_;
    double * Gc_m_;    
    double * Gp_m_x_;
    double * Ge_m_x_;
    double * Gc_m_x_;

    // asymptotic values
    double **ge_bar_last_;
    double **gp_p_last_;
    double **gp_m_last_;
    double **ge_p_last_;
    double **ge_m_last_;
    double **gc_bar_last_;
    double **gc_p_last_;
    double **gc_m_last_;
    double **f_bar_last_;
    double *C_min;
    double *C_max;
     
    double * adv_x_;
    double * diff_x_;
    double * elec_mig_x_; 
    double * Flux_x_n_;

    double table_Psi_wall; // dummy variable     
    double *table_sigma_star;
    double *table_lambda_star;
    double **table_f_bar;
    double **table_ge_bar;
    double **table_gc_bar;
    double **table_gp_m;
    double **table_ge_m;
    double **table_gc_m;
    double **table_gp_p;
    double **table_ge_p;
    double **table_gc_p;

    const char *filename_1_;
    const char *filename_2_;
    const char *filename_3_;
    const char *filename_4_;
    int     sigma_power_;
    int     lambda_power_;
    double  sigma_star_min_;
    double  lambda_star_min_;
    int     max_sigma_star_counter_;
    int     max_lambda_star_counter_;
    int     Num_Channel_;
    int     Num_Reservoir_;
    int     *Channel_Num_Cells_;
    int     Total_Cells_;       //including reservoir cells
    int     Total_Faces_;       //only faces for interial cells
    double  lambda_ref_;
    double  Pe_;
    double  gp_[2];
    double  aa;
    double  bb;
    double  dum1;
    double  dum2;
    int     cell_counter;
    int     face_counter;
    int     channel_counter;
       
    double *Reservoir_Net_Flux_;
    double *Channel_Inlet_Flux_;
    double *Channel_Outlet_Flux_;

    double **three_inlet_flux_;
    double **three_outlet_flux_;

    //Arrays used for implicit solver
    double *implicit_rhs;
    double *implicit_D;
    double *implicit_LD;
    double *implicit_UD;
    double *reservoir_dc;
    double implicit_relaxation;

    //prevent copying and assignment
    Table(const Table &);
    Table & operator= (const Table &);
};
#endif
