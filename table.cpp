#include"table.hpp"

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<stdio.h>

 void Table::setup(int Num_Channel, int Num_Reservoir, int *Channel_Num_Cells, int *Channel_type,
 int *Reservoir_type, double lambda_ref, double Pe, double slit_gp, double circle_gp, 
 const char* filename_1, const char* filename_2,const char* filename_3, const char* filename_4, 
 unsigned max_sigma_star_counter, unsigned max_lambda_star_counter, double sigma_star_min, 
 double lambda_star_min){

  implicit_relaxation = 1.;

  //set some values:
  filename_1_ = filename_1;

  filename_2_ = filename_2;

  filename_3_ = filename_3;

  filename_4_ = filename_4;

  max_sigma_star_counter_ = max_sigma_star_counter;
    
  max_lambda_star_counter_ = max_lambda_star_counter;
   
  sigma_star_min_ = sigma_star_min;
    
  lambda_star_min_ = lambda_star_min;

  lambda_ref_ = lambda_ref;

  Pe_=Pe;
 
  gp_[0] = slit_gp;

  gp_[1] = circle_gp;
    
  Num_Channel_=Num_Channel;
  
  Num_Reservoir_=Num_Reservoir;

  Total_Cells_=0;

  Total_Faces_=0;

  Channel_Num_Cells_=new int[Num_Channel_];
 
  for(int i=0; i<Num_Channel_; i++)
    {
      Channel_Num_Cells_[i]=Channel_Num_Cells[i];

      Total_Cells_+=Channel_Num_Cells[i]+2;

      Total_Faces_+=Channel_Num_Cells[i]+1;
    }
 
  // memory allocation

  Channel_type_ = new int[Num_Channel];
  Reservoir_type_ = new int[Num_Reservoir];

  for(int i=0; i<Num_Channel_; i++)
    {
      Channel_type_[i] = Channel_type[i];
    }

  for(int i=0; i<Num_Reservoir_; i++)
    {
      Reservoir_type_[i] = Reservoir_type[i];
    }
    
  table_sigma_star=new double[max_sigma_star_counter_*max_lambda_star_counter_];
  table_lambda_star=new double[max_sigma_star_counter_*max_lambda_star_counter_];
  table_f_bar=new double*[max_sigma_star_counter_*max_lambda_star_counter_];
  table_ge_bar=new double*[max_sigma_star_counter_*max_lambda_star_counter_];
  table_gc_bar=new double*[max_sigma_star_counter_*max_lambda_star_counter_];
  table_gp_m=new double*[max_sigma_star_counter_*max_lambda_star_counter_];
  table_ge_m=new double*[max_sigma_star_counter_*max_lambda_star_counter_];
  table_gc_m=new double*[max_sigma_star_counter_*max_lambda_star_counter_];
  table_gp_p=new double*[max_sigma_star_counter_*max_lambda_star_counter_];
  table_ge_p=new double*[max_sigma_star_counter_*max_lambda_star_counter_];
  table_gc_p=new double*[max_sigma_star_counter_*max_lambda_star_counter_];

  for(int i=0; i<max_sigma_star_counter_*max_lambda_star_counter_; i++)
    {
      table_f_bar[i] = new double[2];
      table_ge_bar[i] = new double[2];
      table_gc_bar[i] = new double[2];
      table_gp_m[i] = new double[2]; 
      table_ge_m[i] = new double[2];
      table_gc_m[i] = new double[2]; 
      table_gp_p[i] = new double[2];
      table_ge_p[i] = new double[2];
      table_gc_p[i] = new double[2];
    }
 
  // Note:values in cell centers
  f_bar_=new double[Total_Cells_];
  ge_bar_=new double[Total_Cells_];
  gc_bar_=new double[Total_Cells_];
  gp_m_=new double[Total_Cells_];
  ge_m_=new double[Total_Cells_];
  gc_m_=new double[Total_Cells_];
  gp_p_=new double[Total_Cells_];
  ge_p_=new double[Total_Cells_];
  gc_p_=new double[Total_Cells_];

  // extra arrays
  gcbar_over_fbar_=new double[Total_Cells_];
  gcp_over_fbar_=new double[Total_Cells_];
  gcm_over_fbar_=new double[Total_Cells_];
  sqrt_C_over_fbar_=new double[Total_Cells_];
  fbar_sqrt_C_=new double[Total_Cells_];

  Gp_m_=new double[Total_Cells_];
  Ge_m_=new double[Total_Cells_];
  Gc_m_=new double[Total_Cells_];
  
  Gp_m_x_=new double[Total_Faces_];
  Ge_m_x_=new double[Total_Faces_];
  Gc_m_x_=new double[Total_Faces_];

  A1_=new double[Total_Cells_];
  A2_=new double[Total_Cells_];
  B1_=new double[Total_Cells_];
  B2_=new double[Total_Cells_];
  f1_1_=new double[Total_Cells_];
  f2_1_=new double[Total_Cells_];
     
  RHS_=new double[2*Total_Cells_];

  //Note:values in faces
  A1_x_=new double[Total_Faces_];
  A2_x_=new double[Total_Faces_];
  B1_x_=new double[Total_Faces_];
  B2_x_=new double[Total_Faces_];
  f1_x_=new double[Total_Faces_];
  f2_x_=new double[Total_Faces_];
  sqrt_Cf_x_=new double[Total_Faces_];
  dCfdx=new double[Total_Faces_];

  adv_x_=new double[Total_Faces_];
  diff_x_=new double[Total_Faces_];
  elec_mig_x_=new double[Total_Faces_];
  Flux_x_n_=new double[Total_Faces_];
  S_x_=new double[Total_Faces_];
    
  Reservoir_Net_Flux_=new double[Num_Reservoir_];
  Channel_Inlet_Flux_=new double [Num_Channel_];
  Channel_Outlet_Flux_=new double [Num_Channel_];

  three_inlet_flux_ = new double* [Num_Channel_];
  three_outlet_flux_ = new double* [Num_Channel_];
  for(int i=0; i<Num_Channel_; i++)
    {
      three_inlet_flux_[i] = new double [3];
      three_outlet_flux_[i] = new double [3];
    }

  //implicit solver arrays:
  implicit_rhs = new double[Total_Cells_*3];
  implicit_D   = new double[Total_Cells_];
  implicit_LD  = new double[Total_Cells_-1];
  implicit_UD  = new double[Total_Cells_-1];

  reservoir_dc = new double[Num_Reservoir_];  

  // asymptotic values
  ge_bar_last_ = new double*[max_sigma_star_counter_];
  gp_p_last_ = new double*[max_sigma_star_counter_];
  gp_m_last_ = new double*[max_sigma_star_counter_];
  ge_p_last_ = new double*[max_sigma_star_counter_];
  ge_m_last_ = new double*[max_sigma_star_counter_];
  gc_bar_last_ = new double*[max_sigma_star_counter_];
  gc_p_last_ = new double*[max_sigma_star_counter_];
  gc_m_last_ = new double*[max_sigma_star_counter_];
  f_bar_last_ = new double*[max_sigma_star_counter_];
  C_min = new double[Total_Cells_];
  C_max = new double[Total_Cells_];

  for(int i=0; i<max_sigma_star_counter_; i++)
    {

      ge_bar_last_[i] = new double[2];
      gp_p_last_[i] = new double[2];
      gp_m_last_[i] = new double[2];
      ge_p_last_[i] = new double[2];
      ge_m_last_[i] = new double[2];
      gc_bar_last_[i] = new double[2];
      gc_p_last_[i] = new double[2];
      gc_m_last_[i] = new double[2];
      f_bar_last_[i] = new double[2];
    }

  }

Table::~Table(){
    
  ////Deleting Memory Allocation
  delete [] Channel_type_;
  delete [] Reservoir_type_;

  for(int i=0 ; i<max_sigma_star_counter_*max_lambda_star_counter_; i++) {
  
    delete [] table_f_bar[i];
    delete [] table_ge_bar[i];
    delete [] table_gc_bar[i];
    delete [] table_gp_m[i];
    delete [] table_ge_m[i];
    delete [] table_gc_m[i];
    delete [] table_gp_p[i];
    delete [] table_ge_p[i];
    delete [] table_gc_p[i];
  }

  delete [] table_sigma_star;
  delete [] table_lambda_star;
  delete [] table_f_bar;
  delete [] table_ge_bar;
  delete [] table_gc_bar;
  delete [] table_gp_m;
  delete [] table_ge_m;
  delete [] table_gc_m;
  delete [] table_gp_p;
  delete [] table_ge_p;
  delete [] table_gc_p;

  delete [] f_bar_;
  delete [] ge_bar_;
  delete [] gc_bar_;
  delete [] gp_m_;
  delete [] ge_m_;
  delete [] gc_m_;
  delete [] gp_p_;
  delete [] ge_p_;
  delete [] gc_p_;

  delete [] Gp_m_;
  delete [] Ge_m_;
  delete [] Gc_m_;
  delete [] Gp_m_x_;
  delete [] Ge_m_x_;
  delete [] Gc_m_x_;

  delete [] A1_;
  delete [] B1_;
  delete [] A2_;
  delete [] B2_;
  delete [] f1_1_;
  delete [] f2_1_;
  delete [] RHS_;
  
  delete [] A1_x_;
  delete [] A2_x_;
  delete [] B1_x_;
  delete [] B2_x_;
  delete [] f1_x_;
  delete [] f2_x_;

  delete [] adv_x_;
  delete [] diff_x_;
  delete [] elec_mig_x_;
  delete [] Flux_x_n_;
  delete [] S_x_;
  
  delete [] Reservoir_Net_Flux_;
  delete [] Channel_Num_Cells_;

  delete [] Channel_Inlet_Flux_;
  delete [] Channel_Outlet_Flux_;

  for(int i=0; i<Num_Channel_; i++)
    {
      delete [] three_inlet_flux_[i];
      delete [] three_outlet_flux_[i];
    }

  delete [] three_inlet_flux_;
  delete [] three_outlet_flux_;

  // implicit solver arrays:
  delete [] implicit_rhs;
  delete [] implicit_D;
  delete [] implicit_LD;
  delete [] implicit_UD;
  delete [] reservoir_dc;

  // asympotic values
  for(int i=0; i<max_sigma_star_counter_; i++)
    {

      delete [] ge_bar_last_[i];
      delete [] gp_p_last_[i];
      delete [] gp_m_last_[i];
      delete [] ge_p_last_[i];
      delete [] ge_m_last_[i];
      delete [] gc_bar_last_[i];
      delete [] gc_p_last_[i];
      delete [] gc_m_last_[i];
      delete [] f_bar_last_[i];
    }
  delete [] ge_bar_last_;
  delete [] gp_p_last_;
  delete [] gp_m_last_;
  delete [] ge_p_last_;
  delete [] ge_m_last_;
  delete [] f_bar_last_;
  delete [] gc_bar_last_;
  delete [] gc_p_last_;
  delete [] gc_m_last_;
  delete [] C_min; 
  delete [] C_max;

  // extra arrays
  delete [] gcbar_over_fbar_;
  delete [] gcp_over_fbar_;
  delete [] gcm_over_fbar_;
  delete [] sqrt_C_over_fbar_;
  delete [] fbar_sqrt_C_;

  delete [] sqrt_Cf_x_;
  delete [] dCfdx;
}

void Table::Compute_asymptotic(double *lambda_for_table)
{
  for(int i=0; i<max_sigma_star_counter_; i++)
    {
      for(int j=0; j<2; j++)
	{
	  ge_bar_last_[i][j] = table_ge_bar[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j];
	  gp_p_last_[i][j] = table_gp_p[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j];
	  gp_m_last_[i][j] = table_gp_m[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j];
	  ge_p_last_[i][j] = table_ge_p[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j];
	  ge_m_last_[i][j] = table_ge_m[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j]; 
	  f_bar_last_[i][j] = table_f_bar[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j];
	  gc_bar_last_[i][j] = table_gc_bar[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j];
	  gc_p_last_[i][j] = table_gc_p[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j];
	  gc_m_last_[i][j] = table_gc_m[i*max_lambda_star_counter_+max_lambda_star_counter_-1][j];
	}
      
     }

  for(int i=0;i<Total_Cells_;i++)
    {
     C_min[i] = (lambda_for_table[i]/table_lambda_star[max_lambda_star_counter_-1])*(lambda_for_table[i]/table_lambda_star[max_lambda_star_counter_-1]);
     
     C_max[i] = (lambda_for_table[i]/table_lambda_star[0])*(lambda_for_table[i]/table_lambda_star[0]);
  }
  
  return;
}
     
void Table::Import_Table_1(double *lambda_for_table){

  //Importing from input files slit and circle: 
  ifstream importing1(filename_1_);
  ifstream importing3(filename_3_);

    for(int i=0; i<max_sigma_star_counter_*max_lambda_star_counter_;i++)
    {
        importing1>>table_sigma_star[i]>> table_lambda_star[i]>>table_f_bar[i][0]>>

        table_Psi_wall>>table_ge_bar[i][0]>> table_gc_bar[i][0]>>table_gp_m[i][0]>>

	  table_ge_m[i][0]>>table_gc_m[i][0]>>table_gp_p[i][0]>>table_ge_p[i][0]>>table_gc_p[i][0];


        importing3>>table_sigma_star[i]>> table_lambda_star[i]>>table_f_bar[i][1]>>

        table_Psi_wall>>table_ge_bar[i][1]>> table_gc_bar[i][1]>>table_gp_m[i][1]>>

	  table_ge_m[i][1]>>table_gc_m[i][1]>>table_gp_p[i][1]>>table_ge_p[i][1]>>table_gc_p[i][1];

    }
    importing1.close();
    importing3.close();

    Compute_asymptotic(lambda_for_table);

    return;
}
     
void Table::Import_Table_2(double *lambda_for_table){

  //Importing from input files slit and circle: 
  ifstream importing2(filename_2_);
  ifstream importing4(filename_4_);

    for(int i=0; i<max_sigma_star_counter_*max_lambda_star_counter_;i++)
    {
        importing2>>table_sigma_star[i]>> table_lambda_star[i]>>table_f_bar[i][0]>>

        table_Psi_wall>>table_ge_bar[i][0]>> table_gc_bar[i][0]>>table_gp_m[i][0]>>

	table_ge_m[i][0]>>table_gc_m[i][0]>>table_gp_p[i][0]>>table_ge_p[i][0]>>table_gc_p[i][0];


	importing4>>table_sigma_star[i]>> table_lambda_star[i]>>table_f_bar[i][1]>>

        table_Psi_wall>>table_ge_bar[i][1]>> table_gc_bar[i][1]>>table_gp_m[i][1]>>

	table_ge_m[i][1]>>table_gc_m[i][1]>>table_gp_p[i][1]>>table_ge_p[i][1]>>table_gc_p[i][1];
   }
    importing2.close();
    importing4.close();

    Compute_asymptotic(lambda_for_table);
    return;
}
     
void Table::Find_From_Table(double sigma_star, double lambda, double C_bar, int type, int idx){

 double lambda_star;

if(floor(abs(sigma_star)/abs(sigma_star_min_))==0)
{
  //cout<<"NOTE: No surface charge on the surface!"<<endl;

	 f_bar_[idx] = 1.;
	 ge_bar_[idx] = 0.;
	 gc_bar_[idx] = 0;
	 gp_m_[idx] = 0;
	 ge_m_[idx] = 0.;
	 gc_m_[idx] = 0.;
	 gp_p_[idx] = 0.;
	 ge_p_[idx] = 0.;
	 gc_p_[idx] = 0.;

	 gcbar_over_fbar_[idx] = gc_bar_[idx]/(C_bar/f_bar_[idx]);

         gcp_over_fbar_[idx] = gc_p_[idx]/(C_bar/f_bar_[idx]);

         gcm_over_fbar_[idx] = gc_m_[idx]/(C_bar/f_bar_[idx]);

         sqrt_C_over_fbar_[idx] = C_bar/(f_bar_[idx]*f_bar_[idx]);

         fbar_sqrt_C_[idx] = f_bar_[idx]*sqrt(C_bar);
          
}
 else{
	// asymptotic relations
	if(C_bar < C_min[idx])
	  {
	    //cout<<"NOTE: LARGE lambda_star: C["<<i<<"]="<<C_bar[i]<<endl;
     
	    sigma_power_=int(10*log10(sigma_star/sigma_star_min_));
            if(sigma_power_ < 0 ) sigma_power_ = 0;

	    aa=(log10(abs(sigma_star))-log10(abs(table_sigma_star[max_lambda_star_counter_*sigma_power_])))/(log10(abs(table_sigma_star[max_lambda_star_counter_*(sigma_power_+1)]))- log10(abs(table_sigma_star[max_lambda_star_counter_*sigma_power_])));

	    ge_bar_[idx] = aa*(log10(abs(ge_bar_last_[sigma_power_+1][type]))-log10(abs(ge_bar_last_[sigma_power_][type]))) + log10(abs(ge_bar_last_[sigma_power_][type]));
	    ge_bar_[idx] = -pow(10, ge_bar_[idx]);

	    gp_p_[idx] = aa*(log10(abs(gp_p_last_[sigma_power_+1][type]))-log10(abs(gp_p_last_[sigma_power_][type]))) + log10(abs(gp_p_last_[sigma_power_][type]));
	    gp_p_[idx] = pow(10, gp_p_[idx]);

	    ge_p_[idx] = aa*(log10(abs(ge_p_last_[sigma_power_+1][type]))-log10(abs(ge_p_last_[sigma_power_][type]))) + log10(abs(ge_p_last_[sigma_power_][type]));
	    ge_p_[idx] = pow(10, ge_p_[idx]);

	    gp_m_[idx] = aa*(log10(abs(gp_m_last_[sigma_power_+1][type]))-log10(abs(gp_m_last_[sigma_power_][type]))) + log10(abs(gp_m_last_[sigma_power_][type]));
	    gp_m_[idx] = -pow(10, gp_m_[idx]);
	    gp_m_[idx] = gp_m_[idx]*(C_bar/C_min[idx]);
	   
	    ge_m_[idx] = aa*(log10(abs(ge_m_last_[sigma_power_+1][type]))-log10(abs(ge_m_last_[sigma_power_][type]))) + log10(abs(ge_m_last_[sigma_power_][type]));
	    ge_m_[idx] = -pow(10, ge_m_[idx]);
	    ge_m_[idx] = ge_m_[idx]*(C_bar/C_min[idx]);
	   
	    gc_bar_[idx] = aa*(log10(abs(gc_bar_last_[sigma_power_+1][type]))-log10(abs(gc_bar_last_[sigma_power_][type]))) + log10(abs(gc_bar_last_[sigma_power_][type]));
	    gc_bar_[idx] = -pow(10, gc_bar_[idx]);

	     gc_p_[idx] = aa*(log10(abs(gc_p_last_[sigma_power_+1][type]))-log10(abs(gc_p_last_[sigma_power_][type]))) + log10(abs(gc_p_last_[sigma_power_][type]));
	     gc_p_[idx] = pow(10, gc_p_[idx]);

	     gc_m_[idx] = aa*(log10(abs(gc_m_last_[sigma_power_+1][type]))-log10(abs(gc_m_last_[sigma_power_][type]))) + log10(abs(gc_m_last_[sigma_power_][type]));
	     gc_m_[idx] = -pow(10, gc_m_[idx]);

	     f_bar_[idx] = aa*(log10(abs(f_bar_last_[sigma_power_+1][type]))-log10(abs(f_bar_last_[sigma_power_][type]))) + log10(abs(f_bar_last_[sigma_power_][type]));
	     f_bar_[idx] = pow(10, f_bar_[idx]);

	     gcbar_over_fbar_[idx] = gc_bar_[idx]/(C_min[idx]/f_bar_[idx]);

	     gcp_over_fbar_[idx] = gc_p_[idx]/(C_min[idx]/f_bar_[idx]);

	     gcm_over_fbar_[idx] = (gc_m_[idx]/(C_min[idx]/f_bar_[idx]))*(C_bar/C_min[idx]); 

	     sqrt_C_over_fbar_[idx] =  C_min[idx]/(f_bar_[idx]*f_bar_[idx]);

	     fbar_sqrt_C_[idx] = f_bar_[idx]*sqrt(C_min[idx])*(C_bar/C_min[idx]);
	    
	  }
        else if(C_bar > C_max[idx])
        {
	  //cout<<"NOTE: SMALL lambda_star: C["<<i<<"]="<<C_bar[i]<<endl;
          // Extrapolation
    
	    lambda_star = lambda/sqrt(C_bar);
 
	    sigma_power_= int(10*log10(sigma_star/sigma_star_min_));
            if(sigma_power_ < 0 ) sigma_power_ = 0;

	    aa=(log10(abs(sigma_star))-log10(abs(table_sigma_star[max_lambda_star_counter_*sigma_power_])))/(log10(abs(table_sigma_star[max_lambda_star_counter_*(sigma_power_+1)]))- log10(abs(table_sigma_star[max_lambda_star_counter_*sigma_power_])));

	    bb=(log10(table_lambda_star[max_lambda_star_counter_*sigma_power_])-log10(lambda_star))/(log10(table_lambda_star[max_lambda_star_counter_*sigma_power_+1])-log10(table_lambda_star[max_lambda_star_counter_*sigma_power_]));
        
	       
        //f_bar
        dum1=aa*(log10(abs(table_f_bar[max_lambda_star_counter_*(sigma_power_+1)][type]))-
                 
                 log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_][type])))+
        
           log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_][type]));
        
        dum2=aa*(log10(abs(table_f_bar[max_lambda_star_counter_*(sigma_power_+1)+1][type]))
                 
                 -log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+1][type])))+
        
           log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+1][type]));
        
        f_bar_[idx]= dum1 - bb*(dum2-dum1);
        
        f_bar_[idx]=pow(10,f_bar_[idx]);
        
	//ge_bar
        dum1=aa*(log10(abs(table_ge_bar[max_lambda_star_counter_*(sigma_power_+1)][type]))-
                 
                 log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_][type])))+
        
           log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_][type]));
        
        dum2=aa*(log10(abs(table_ge_bar[max_lambda_star_counter_*(sigma_power_+1)+1][type]))-
                 
                 log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+1][type])))+
        
           log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+1][type]));
        
        ge_bar_[idx]=dum1 - bb*(dum2-dum1);
        
        ge_bar_[idx]=-pow(10,ge_bar_[idx]);
        
        //gc_bar
        dum1=aa*(log10(abs(table_gc_bar[max_lambda_star_counter_*(sigma_power_+1)][type]))-
                 
                 log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_][type])))+
        
           log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_][type]));

        dum2=aa*(log10(abs(table_gc_bar[max_lambda_star_counter_*(sigma_power_+1)+1][type]))-
                 
                 log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+1][type])))+
        
           log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+1][type]));
        
        gc_bar_[idx]=dum1 - bb*(dum2-dum1);
        
        gc_bar_[idx]= -pow(10,gc_bar_[idx]);
        
        
        //gp_m
        dum1=aa*(log10(abs(table_gp_m[max_lambda_star_counter_*(sigma_power_+1)][type]))-
                 
                 log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_][type])))+
        
           log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_][type]));
        
        dum2=aa*(log10(abs(table_gp_m[max_lambda_star_counter_*(sigma_power_+1)+1][type]))
                 
                 -log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+1][type])))+
        
           log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+1][type]));
        
        gp_m_[idx]=dum1 - bb*(dum2-dum1);
        
        gp_m_[idx]= -pow(10,gp_m_[idx]);

        //ge_m
        dum1=aa*(log10(abs(table_ge_m[max_lambda_star_counter_*(sigma_power_+1)][type]))-
                 
                 log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_][type])))+
        
           log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_][type]));
        
        dum2=aa*(log10(abs(table_ge_m[max_lambda_star_counter_*(sigma_power_+1)+1][type]))-
                 
                 log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+1][type])))+
        
           log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+1][type]));
        
        ge_m_[idx]=dum1 - bb*(dum2-dum1);
        
        ge_m_[idx]= -pow(10,ge_m_[idx]);
        
        
        //gc_m
        dum1=aa*(log10(abs(table_gc_m[max_lambda_star_counter_*(sigma_power_+1)][type]))-
                 
                 log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_][type])))+
        
           log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_][type]));

        dum2=aa*(log10(abs(table_gc_m[max_lambda_star_counter_*(sigma_power_+1)+1][type]))-
                 
                 log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+1][type])))+
        
           log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+1][type]));
        
        gc_m_[idx]=dum1 - bb*(dum2-dum1);
        
        gc_m_[idx]= -pow(10,gc_m_[idx]);
        
        //gp_p
        dum1=aa*(log10(abs(table_gp_p[max_lambda_star_counter_*(sigma_power_+1)][type]))
                 
                 -log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_][type])))+
        
           log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_][type]));
        
        dum2=aa*(log10(abs(table_gp_p[max_lambda_star_counter_*(sigma_power_+1)+1][type]))-
                 
                 log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+1][type])))+
        
           log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+1][type]));
        
        gp_p_[idx]=dum1 - bb*(dum2-dum1);
        
        gp_p_[idx]= pow(10,gp_p_[idx]);

        //ge_p
        dum1=aa*(log10(abs(table_ge_p[max_lambda_star_counter_*(sigma_power_+1)][type]))-
                 
                 log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_][type])))+
        
	  log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_][type]));
        
        dum2=aa*(log10(abs(table_ge_p[max_lambda_star_counter_*(sigma_power_+1)+1][type]))-
                 
                 log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+1][type])))+
        
	  log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+1][type]));
        
        ge_p_[idx]=dum1 - bb*(dum2-dum1);
        
        ge_p_[idx]= pow(10,ge_p_[idx]);
        
        //gc_p
        dum1=aa*(log10(abs(table_gc_p[max_lambda_star_counter_*(sigma_power_+1)][type]))-
                 
                 log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_][type])))+
        
	  log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_][type]));

         dum2=aa*(log10(abs(table_gc_p[max_lambda_star_counter_*(sigma_power_+1)+1][type]))-
                 
                 log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+1][type])))+
        
	   log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+1][type]));
        
	 gc_p_[idx]=dum1 - bb*(dum2-dum1);
        
	 gc_p_[idx]= pow(10,gc_p_[idx]);


	 gcbar_over_fbar_[idx] = gc_bar_[idx]/(C_bar/f_bar_[idx]);

	 gcp_over_fbar_[idx] = gc_p_[idx]/(C_bar/f_bar_[idx]);

	 gcm_over_fbar_[idx] = gc_m_[idx]/(C_bar/f_bar_[idx]);

	 sqrt_C_over_fbar_[idx] =  C_bar/(f_bar_[idx]*f_bar_[idx]);

	 fbar_sqrt_C_[idx] = f_bar_[idx]*sqrt(C_bar);
	    
	  }
		
	else
	  {
	    
        lambda_star = lambda/sqrt(C_bar);

        sigma_power_ = int(10*log10(sigma_star/sigma_star_min_));
        
        lambda_power_ = int(10*log10(lambda_star/lambda_star_min_));

	if(sigma_power_ < 0 ) sigma_power_ = 0;

	   aa=(log10(abs(sigma_star))-log10(abs(table_sigma_star[max_lambda_star_counter_*sigma_power_+lambda_power_])))/(log10(abs(table_sigma_star[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))- log10(abs(table_sigma_star[max_lambda_star_counter_*sigma_power_+lambda_power_])));

	    bb=(log10(lambda_star)-log10(table_lambda_star[max_lambda_star_counter_*sigma_power_+lambda_power_]))/(log10(table_lambda_star[max_lambda_star_counter_*sigma_power_+lambda_power_+1])-log10(table_lambda_star[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
	       
        //f_bar
        dum1=aa*(log10(abs(table_f_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))-
                 
                 log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
           log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));
        
        dum2=aa*(log10(abs(table_f_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))
                 
                 -log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
           log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
        f_bar_[idx]=bb*(dum2-dum1)+dum1;
        
        f_bar_[idx]=pow(10,f_bar_[idx]);
        
	//ge_bar
        dum1=aa*(log10(abs(table_ge_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))-
                 
                 log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
           log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));
        
        dum2=aa*(log10(abs(table_ge_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))-
                 
                 log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
           log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
        ge_bar_[idx]=bb*(dum2-dum1)+dum1;
        
        ge_bar_[idx]=-pow(10,ge_bar_[idx]);
        
        //gc_bar
        dum1=aa*(log10(abs(table_gc_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))-
                 
                 log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
           log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));

        dum2=aa*(log10(abs(table_gc_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))-
                 
                 log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
           log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
        gc_bar_[idx]=bb*(dum2-dum1)+dum1;
        
        gc_bar_[idx]= -pow(10,gc_bar_[idx]);
        
        
        //gp_m
        dum1=aa*(log10(abs(table_gp_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))-
                 
                 log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
           log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));
        
        dum2=aa*(log10(abs(table_gp_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))
                 
                 -log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
           log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
        gp_m_[idx]=bb*(dum2-dum1)+dum1;
        
        gp_m_[idx]= -pow(10,gp_m_[idx]);

        //ge_m
        dum1=aa*(log10(abs(table_ge_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))-
                 
                 log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
           log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));
        
        dum2=aa*(log10(abs(table_ge_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))-
                 
                 log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
           log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
        ge_m_[idx]=bb*(dum2-dum1)+dum1;
        
        ge_m_[idx]= -pow(10,ge_m_[idx]);
        
        
        //gc_m
        dum1=aa*(log10(abs(table_gc_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))-
                 
                 log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
           log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));

        dum2=aa*(log10(abs(table_gc_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))-
                 
                 log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
           log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
        gc_m_[idx]=bb*(dum2-dum1)+dum1;
        
        gc_m_[idx]= -pow(10,gc_m_[idx]);
        
        //gp_p
        dum1=aa*(log10(abs(table_gp_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))
                 
                 -log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
           log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));
        
        dum2=aa*(log10(abs(table_gp_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))-
                 
                 log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
           log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
        gp_p_[idx]=bb*(dum2-dum1)+dum1;
        
        gp_p_[idx]= pow(10,gp_p_[idx]);

        //ge_p
        dum1=aa*(log10(abs(table_ge_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))-
                 
                 log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
	  log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));
        
        dum2=aa*(log10(abs(table_ge_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))-
                 
                 log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
	  log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
        ge_p_[idx]=bb*(dum2-dum1)+dum1;
        
        ge_p_[idx]= pow(10,ge_p_[idx]);
        
        //gc_p
        dum1=aa*(log10(abs(table_gc_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_][type]))-
                 
                 log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+lambda_power_][type])))+
        
	  log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+lambda_power_][type]));

         dum2=aa*(log10(abs(table_gc_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1][type]))-
                 
                 log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type])))+
        
	   log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1][type]));
        
	 gc_p_[idx]=bb*(dum2-dum1)+dum1;
        
	 gc_p_[idx]= pow(10,gc_p_[idx]);


	 gcbar_over_fbar_[idx] = gc_bar_[idx]/(C_bar/f_bar_[idx]);

	 gcp_over_fbar_[idx] = gc_p_[idx]/(C_bar/f_bar_[idx]);

	 gcm_over_fbar_[idx] = gc_m_[idx]/(C_bar/f_bar_[idx]);

	 sqrt_C_over_fbar_[idx] = C_bar/(f_bar_[idx]*f_bar_[idx]);

	 fbar_sqrt_C_[idx] = f_bar_[idx]*sqrt(C_bar);

       } // End of else for C<0
    
 }//End of Else

return;
}

void Table::Find_Initial_C_bar(double *sigma_star, double *lambda_for_table, double *C0, int **Channel_End_Reservoir, double *C_bar, double *Reservoir_C_bar, double time){
    
    //Importing table based on sigma_star and lambda_0
    Import_Table_1(lambda_for_table);

    cell_counter = 0;
   
    for(int i = 0; i<Num_Channel_; i++){
      
      Find_From_Table(sigma_star[cell_counter], lambda_for_table[cell_counter], C0[cell_counter], Reservoir_type_[Channel_End_Reservoir[i][0]], cell_counter);

      for(int j = 1; j <= Channel_Num_Cells_[i]; j++){
	
	Find_From_Table(sigma_star[cell_counter+j], lambda_for_table[cell_counter+j], C0[cell_counter+j],Channel_type_[i], cell_counter+j);
      }

      cell_counter += Channel_Num_Cells_[i]+2;

      Find_From_Table(sigma_star[cell_counter-1], lambda_for_table[cell_counter-1], C0[cell_counter-1], Reservoir_type_[Channel_End_Reservoir[i][1]], cell_counter-1);
   
    }

    //finding initial C_bar
    for(int i=0;i<Total_Cells_;i++) {
	C_bar[i] = C0[i]*f_bar_[i];
    }
   
    //Storing C_bar of reservoirs in a seperate array used for visualization    
    cell_counter=0;
    for(int i=0; i<Num_Channel_; i++)
    {
        Reservoir_C_bar[Channel_End_Reservoir[i][0]]=C_bar[cell_counter];
         
        cell_counter+=Channel_Num_Cells_[i]+2;
        
        Reservoir_C_bar[Channel_End_Reservoir[i][1]]=C_bar[cell_counter-1];
    }
  
    //Removing first table data and importing new table
    //based on sigma_star and lambda_star, this table includes dlog_f_bar!
    Import_Table_2(lambda_for_table);

    return;
}

void Table::calculate_A_B_f_RHS(double *C_bar, double *C0, double *Cs, double *lambda, double *nondimensional_Sp, double *dx, double *Art_Diff, int **Channel_End_Reservoir){
    
  ////in cell centers
  int type = 0;
  cell_counter = 0;

  for(int i=0; i<Num_Channel_; i++){

    type = Reservoir_type_[Channel_End_Reservoir[i][0]];
    
    A1_[cell_counter] = Pe_/(2*lambda[cell_counter]*lambda[cell_counter])*gp_[type];
        
    B1_[cell_counter] = Pe_*ge_bar_[cell_counter];
        
    f1_1_[cell_counter] = -Pe_/(2*lambda[cell_counter]*lambda[cell_counter])*0.5*gcbar_over_fbar_[cell_counter];

    A2_[cell_counter] = Pe_/(2*lambda[cell_counter]*lambda[cell_counter])*Cs[cell_counter]*(gp_[type]+gp_p_[cell_counter]-gp_m_[cell_counter]);
	        
    B2_[cell_counter] = -(2*C_bar[cell_counter]+Cs[cell_counter])+ Pe_*Cs[cell_counter]*(ge_bar_[cell_counter]+ge_p_[cell_counter]-ge_m_[cell_counter]);
        
    f2_1_[cell_counter] = -Pe_/(2*lambda[cell_counter]*lambda[cell_counter])*Cs[cell_counter]*0.5*(gcbar_over_fbar_[cell_counter]+gcp_over_fbar_[cell_counter]-gcm_over_fbar_[cell_counter]) - 1./sqrt_C_over_fbar_[cell_counter];

    Gp_m_[cell_counter] = Pe_/(2*lambda[cell_counter]*lambda[cell_counter])*gp_m_[cell_counter]*Cs[cell_counter];
        
    Ge_m_[cell_counter] = Pe_*ge_m_[cell_counter]*Cs[cell_counter];

    Gc_m_[cell_counter] = Pe_/(2*lambda[cell_counter]*lambda[cell_counter])*0.5*gcm_over_fbar_[cell_counter]*Cs[cell_counter];

    type = Channel_type_[i];

    for(int j=1; j<=Channel_Num_Cells_[i]; j++){

        A1_[cell_counter+j] = Pe_/(2*lambda[cell_counter+j]*lambda[cell_counter+j])*gp_[type];
        
        B1_[cell_counter+j] = Pe_*ge_bar_[cell_counter+j];
        
        f1_1_[cell_counter+j] = -Pe_/(2*lambda[cell_counter+j]*lambda[cell_counter+j])*0.5*gcbar_over_fbar_[cell_counter+j];

        A2_[cell_counter+j] = Pe_/(2*lambda[cell_counter+j]*lambda[cell_counter+j])*Cs[cell_counter+j]*(gp_[type]+gp_p_[cell_counter+j]-gp_m_[cell_counter+j]);
	        
        B2_[cell_counter+j] = -(2*C_bar[cell_counter+j]+Cs[cell_counter+j])+ Pe_*Cs[cell_counter+j]*(ge_bar_[cell_counter+j]+ge_p_[cell_counter+j]-ge_m_[cell_counter+j]);
        
        f2_1_[cell_counter+j] = -Pe_/(2*lambda[cell_counter+j]*lambda[cell_counter+j])*Cs[cell_counter+j]*0.5*(gcbar_over_fbar_[cell_counter+j]+gcp_over_fbar_[cell_counter+j]-gcm_over_fbar_[cell_counter+j]) - 1./sqrt_C_over_fbar_[cell_counter+j];

        Gp_m_[cell_counter+j] = Pe_/(2*lambda[cell_counter+j]*lambda[cell_counter+j])*gp_m_[cell_counter+j]*Cs[cell_counter+j];
        
        Ge_m_[cell_counter+j] = Pe_*ge_m_[cell_counter+j]*Cs[cell_counter+j];

	Gc_m_[cell_counter+j] = Pe_/(2*lambda[cell_counter+j]*lambda[cell_counter+j])*0.5*gcm_over_fbar_[cell_counter+j]*Cs[cell_counter+j];
        
     }

    type = Reservoir_type_[Channel_End_Reservoir[i][1]];

    cell_counter += Channel_Num_Cells_[i]+2;

    A1_[cell_counter-1] = Pe_/(2*lambda[cell_counter-1]*lambda[cell_counter-1])*gp_[type];
        
    B1_[cell_counter-1] = Pe_*ge_bar_[cell_counter-1];
        
    f1_1_[cell_counter-1] = -Pe_/(2*lambda[cell_counter-1]*lambda[cell_counter-1])*0.5*gcbar_over_fbar_[cell_counter-1];

    A2_[cell_counter-1] = Pe_/(2*lambda[cell_counter-1]*lambda[cell_counter-1])*Cs[cell_counter-1]*(gp_[type]+gp_p_[cell_counter-1]-gp_m_[cell_counter-1]);
	        
    B2_[cell_counter-1] = -(2*C_bar[cell_counter-1]+Cs[cell_counter-1])+ Pe_*Cs[cell_counter-1]*(ge_bar_[cell_counter-1]+ge_p_[cell_counter-1]-ge_m_[cell_counter-1]);
        
    f2_1_[cell_counter-1] = -Pe_/(2*lambda[cell_counter-1]*lambda[cell_counter-1])*Cs[cell_counter-1]*0.5*(gcbar_over_fbar_[cell_counter-1]+gcp_over_fbar_[cell_counter-1]-gcm_over_fbar_[cell_counter-1]) - 1./sqrt_C_over_fbar_[cell_counter-1];

    Gp_m_[cell_counter-1] = Pe_/(2*lambda[cell_counter-1]*lambda[cell_counter-1])*gp_m_[cell_counter-1]*Cs[cell_counter-1];
        
    Ge_m_[cell_counter-1] = Pe_*ge_m_[cell_counter-1]*Cs[cell_counter-1];

    Gc_m_[cell_counter-1] = Pe_/(2*lambda[cell_counter-1]*lambda[cell_counter-1])*0.5*gcm_over_fbar_[cell_counter-1]*Cs[cell_counter-1];
  }

 
    cell_counter=0; 
    face_counter=0;
    
    ////on faces
    for(int i=0; i<Num_Channel_; i++){
    
        for(int j=0; j<Channel_Num_Cells_[i]+1; j++){

	    	sqrt_Cf_x_[face_counter] = (sqrt_C_over_fbar_[cell_counter+j+1]+sqrt_C_over_fbar_[cell_counter+j])/2.;

	    	diff_x_[face_counter] = (sqrt_C_over_fbar_[cell_counter+j+1]*C_bar[cell_counter+j+1]-sqrt_C_over_fbar_[cell_counter+j]*C_bar[cell_counter+j]) / dx[i];
	    	diff_x_[face_counter] *= Art_Diff[face_counter]; 

	   	 	S_x_[face_counter] = (nondimensional_Sp[cell_counter+j]+nondimensional_Sp[cell_counter+j+1])/2.;

	    	Gp_m_x_[face_counter] = (Gp_m_[cell_counter+j]+Gp_m_[cell_counter+j+1])/2.;

	    	Ge_m_x_[face_counter] = (Ge_m_[cell_counter+j]+Ge_m_[cell_counter+j+1])/2.;

	    	Gc_m_x_[face_counter] = (Gc_m_[cell_counter+j]+Gc_m_[cell_counter+j+1])/2.;
	
            A1_x_[face_counter]=(A1_[cell_counter+j]+A1_[cell_counter+j+1])/2.;
 
            B1_x_[face_counter]=(B1_[cell_counter+j]+B1_[cell_counter+j+1])/2.;

            A2_x_[face_counter]=(A2_[cell_counter+j]+A2_[cell_counter+j+1])/2.;
	           
            B2_x_[face_counter]=(B2_[cell_counter+j]+B2_[cell_counter+j+1])/2.;

	    	// f1 calculation
	    	f1_x_[face_counter] = (f1_1_[cell_counter+j+1]+f1_1_[cell_counter+j])/2.0 * diff_x_[face_counter];
     
            // f2 calculation
	    	f2_x_[face_counter] = (f2_1_[cell_counter+j+1]+f2_1_[cell_counter+j])/2.0 * diff_x_[face_counter];

            face_counter++;
        }

        cell_counter+=Channel_Num_Cells_[i]+2;
    }
    
    ////RHS calculation for ONLY INTERIOR CELLS
    cell_counter=0;
    face_counter=0;
    
    for(int i=0; i< Num_Channel_; i++)
    {
        //setting left boundary element rhs 0, correct values assigned later in poisson solver
        RHS_[cell_counter*2]=0.;
        RHS_[cell_counter*2+1]=0.;
        
        for(int j=1; j<Channel_Num_Cells_[i]+1; j++)
        {
            RHS_[(cell_counter+j)*2]=(S_x_[face_counter+j]*f1_x_[face_counter+j]-S_x_[face_counter+j-1]*f1_x_[face_counter+j-1])/dx[i];
            
            RHS_[(cell_counter+j)*2+1]=(S_x_[face_counter+j]*f2_x_[face_counter+j]-S_x_[face_counter+j-1]*f2_x_[face_counter+j-1])/dx[i];
        }
        
        cell_counter+=Channel_Num_Cells_[i]+2;
        face_counter+=Channel_Num_Cells_[i]+1;
        
        //setting right boundary element rhs 0, correct values assigned later in poisson solver
        RHS_[(cell_counter-1)*2]=0.;
        RHS_[(cell_counter-1)*2+1]=0.;
    }
    
        return;
}

void Table::do_table_reading(double *sigma_star, double *lambda_for_table, double *lambda, double *C_bar, double *C0, double *Reservoir_C0, int **Channel_End_Reservoir, double *Cs, double *nondimensional_Sp, double *dx, double *Art_Diff, double time){

  // Read from table
  cell_counter = 0;
   
  for(int i = 0; i<Num_Channel_; i++){
      
      Find_From_Table(sigma_star[cell_counter], lambda_for_table[cell_counter], C_bar[cell_counter], Reservoir_type_[Channel_End_Reservoir[i][0]], cell_counter);

      for(int j = 1; j <= Channel_Num_Cells_[i]; j++){
    
	Find_From_Table(sigma_star[cell_counter+j], lambda_for_table[cell_counter+j], C_bar[cell_counter+j], Channel_type_[i], cell_counter+j);
      }

      cell_counter += Channel_Num_Cells_[i]+2;

      Find_From_Table(sigma_star[cell_counter-1], lambda_for_table[cell_counter-1], C_bar[cell_counter-1], Reservoir_type_[Channel_End_Reservoir[i][1]], cell_counter-1);

    }

  for(int i=0; i<Total_Cells_; i++) {
    
    C0[i] = C_bar[i]/f_bar_[i];
}
   
  ////NOTE: update Reservoir_C0 with new C0 quantities
  cell_counter = 0;
  for(int i=0; i<Num_Channel_; i++)
    {
        Reservoir_C0[Channel_End_Reservoir[i][0]]=C0[cell_counter];
         
        cell_counter+=Channel_Num_Cells_[i]+2;
        
        Reservoir_C0[Channel_End_Reservoir[i][1]]=C0[cell_counter-1];
    }

  calculate_A_B_f_RHS(C_bar, C0, Cs, lambda, nondimensional_Sp, dx, Art_Diff, Channel_End_Reservoir);
  return;
}

void Table::find_flux_n(double *C_bar, double *C0, double *dPdx, double *dMudx,double *nondimensional_Sp, double *dx, double *U, double time){
   
    ////Velocity and flux in time step n
    face_counter=0;
    cell_counter=0;

    for(int i=0; i< Num_Channel_; i++){
    
      for(int j=0; j<Channel_Num_Cells_[i]+1; j++){
      
      	U[face_counter+j]=A1_x_[face_counter+j]*dPdx[face_counter+j] + B1_x_[face_counter+j]*dMudx[face_counter+j] - f1_x_[face_counter+j];
	
	  	adv_x_[face_counter+j]=U[face_counter+j]*(C_bar[cell_counter+j+1]+C_bar[cell_counter+j])/2. + Gp_m_x_[face_counter+j]* dPdx[face_counter+j] + Ge_m_x_[face_counter+j]*dMudx[face_counter+j] + Gc_m_x_[face_counter+j]*diff_x_[face_counter+j];
          
	  	elec_mig_x_[face_counter+j]=(C_bar[cell_counter+j+1]+C_bar[cell_counter+j])/2.*dMudx[face_counter+j];
            
	  	Flux_x_n_[face_counter+j] = adv_x_[face_counter+j] - 1./sqrt_Cf_x_[face_counter+j]* diff_x_[face_counter+j] + elec_mig_x_[face_counter+j];
            
	  	Flux_x_n_[face_counter+j] = Flux_x_n_[face_counter+j] * S_x_[face_counter+j];
 
        //Storing inlet and outlet fluxes
	    if(j==0){
	    
			Channel_Inlet_Flux_[i]=Flux_x_n_[face_counter+j];
	    }
        else if(j==Channel_Num_Cells_[i]){
        
			Channel_Outlet_Flux_[i]=Flux_x_n_[face_counter+j];
	    }
      }
        
        cell_counter+=Channel_Num_Cells_[i]+2;
        face_counter+=Channel_Num_Cells_[i]+1;
    }

    return;
}

void Table::solve_for_dc(Lapack  &solver, double *C_bar, double *dMudx, double *U, double *nondimensional_Sp, bool *Reservoir_Pressure_Type,  bool *Reservoir_Potential_Type, double *Reservoir_Volume, int **Connectivity, int **Connecting_Channel, int **Channel_End_Reservoir, double *dx, double dt, double *dC_bar)
{
    ////setting everything to zero
    for(int i=0; i<Total_Cells_-1; i++)
    {
      implicit_D[i] = 0.;
      implicit_LD[i] = 0.;
      implicit_UD[i] = 0.;

    }
    implicit_D[Total_Cells_-1] = 0.;

    ////calculating RHS for implicit system
    cell_counter = 0;
    face_counter = 0;

    for(int i=0; i< Num_Channel_; i++)
    {
      implicit_rhs[cell_counter] = 0.;
      implicit_rhs[Total_Cells_+cell_counter] = 1.;
      implicit_rhs[2*Total_Cells_+cell_counter] = 0.;

        for(int j=1; j< Channel_Num_Cells_[i]+1; j++)
        {	             
            implicit_rhs[cell_counter+j]=-(Flux_x_n_[face_counter+j]-Flux_x_n_[face_counter+j-1])/dx[i];

	    implicit_rhs[Total_Cells_+cell_counter+j] = 0.;

	    implicit_rhs[2*Total_Cells_+cell_counter+j] = 0.;       
        }

        cell_counter+=Channel_Num_Cells_[i]+2;
        face_counter+=Channel_Num_Cells_[i]+1;

	implicit_rhs[cell_counter-1] = 0.;
	implicit_rhs[Total_Cells_+cell_counter-1] = 0.;
	implicit_rhs[2*Total_Cells_+cell_counter-1] = 1.;
    }

    ////filling diagonal arrays
    cell_counter = 0;
    face_counter = 0;

    for(int i=0; i<Num_Channel_; i++){

      ////filling main diagonal
      implicit_D[cell_counter] = 1.;

      for(int j=1; j< Channel_Num_Cells_[i]+1; j++){

        //// main diagonal
		implicit_D[cell_counter+j] = implicit_relaxation*S_x_[face_counter+j]/dx[i]*
		(U[face_counter+j]/2.0 + 1./dx[i]*sqrt_C_over_fbar_[cell_counter+j]/sqrt_Cf_x_[face_counter+j] + 1.0/2.0*dMudx[face_counter+j]);

		implicit_D[cell_counter+j] += implicit_relaxation*S_x_[face_counter+j-1]/dx[i]*
		(-U[face_counter+j-1]/2.0 + 1.0/dx[i]*sqrt_C_over_fbar_[cell_counter+j]/sqrt_Cf_x_[face_counter+j-1] - 1.0/2.0*dMudx[face_counter+j-1]);

		implicit_D[cell_counter+j] += nondimensional_Sp[cell_counter+j]/dt;

        //// lower-diagonal
		implicit_LD[cell_counter+j-1] = implicit_relaxation*S_x_[face_counter+j-1]/dx[i]*
		(-U[face_counter+j-1]/2.0 - 1.0/dx[i]*sqrt_C_over_fbar_[cell_counter+j-1]/sqrt_Cf_x_[face_counter+j-1] - 1.0/2.0*dMudx[face_counter+j-1]);

		//// upper-diagonal
		implicit_UD[cell_counter+j] = implicit_relaxation*S_x_[face_counter+j]/dx[i]*
		(U[face_counter+j]/2.0 - 1.0/dx[i]*sqrt_C_over_fbar_[cell_counter+j+1]/sqrt_Cf_x_[face_counter+j]+ 1.0/2.0*dMudx[face_counter+j]);
    }

      	cell_counter += Channel_Num_Cells_[i]+2;
      	face_counter += Channel_Num_Cells_[i]+1;
    
      	implicit_D[cell_counter-1] = 1.;
    }
    
    //// now solving the tri-diagonal systems
    solver.dgtsv(implicit_LD, implicit_D, implicit_UD, implicit_rhs, Total_Cells_ , 3);

    //// find F1, F2 and F3 in the inlet and outlet of the channels
    cell_counter = 0;
    face_counter = 0;
    for(int i=0; i<Num_Channel_; i++){
    
		double inlet_coef = (U[face_counter] + dMudx[face_counter]) ;

		three_inlet_flux_[i][0] = inlet_coef * (implicit_rhs[cell_counter+1]+implicit_rhs[cell_counter])/2.0 - 1.0/sqrt_Cf_x_[face_counter] * (sqrt_C_over_fbar_[cell_counter+1]*implicit_rhs[cell_counter+1] - sqrt_C_over_fbar_[cell_counter]*implicit_rhs[cell_counter])/dx[i];

		three_inlet_flux_[i][0] *= S_x_[face_counter];

		three_inlet_flux_[i][1] = inlet_coef * (implicit_rhs[Total_Cells_+cell_counter+1]+implicit_rhs[Total_Cells_+cell_counter])/2.0 - 1.0/sqrt_Cf_x_[face_counter] * (sqrt_C_over_fbar_[cell_counter+1]*implicit_rhs[Total_Cells_+cell_counter+1] - sqrt_C_over_fbar_[cell_counter]*implicit_rhs[Total_Cells_+cell_counter])/dx[i];
	
		three_inlet_flux_[i][1] *= implicit_relaxation*S_x_[face_counter];

		three_inlet_flux_[i][2] = inlet_coef * (implicit_rhs[2*Total_Cells_+cell_counter+1]+implicit_rhs[2*Total_Cells_+cell_counter])/2.0 - 1.0/sqrt_Cf_x_[face_counter] * (sqrt_C_over_fbar_[cell_counter+1]*implicit_rhs[2*Total_Cells_+cell_counter+1] - sqrt_C_over_fbar_[cell_counter]*implicit_rhs[2*Total_Cells_+cell_counter])/dx[i];

		three_inlet_flux_[i][2] *= implicit_relaxation*S_x_[face_counter];

		cell_counter += Channel_Num_Cells_[i]+2;
		face_counter += Channel_Num_Cells_[i]+1;

		double outlet_coef = (U[face_counter-1] + dMudx[face_counter-1]);

		three_outlet_flux_[i][0] = outlet_coef * (implicit_rhs[cell_counter-1]+implicit_rhs[cell_counter-2])/2.0 - 1.0/sqrt_Cf_x_[face_counter-1] * (sqrt_C_over_fbar_[cell_counter-1]*implicit_rhs[cell_counter-1] - sqrt_C_over_fbar_[cell_counter-2]*implicit_rhs[cell_counter-2])/dx[i] ;

		three_outlet_flux_[i][0] *= S_x_[face_counter-1];

		three_outlet_flux_[i][1] = outlet_coef * (implicit_rhs[Total_Cells_+cell_counter-1]+implicit_rhs[Total_Cells_+cell_counter-2])/2.0 - 1.0/sqrt_Cf_x_[face_counter-1] * (sqrt_C_over_fbar_[cell_counter-1]*implicit_rhs[Total_Cells_+cell_counter-1] - sqrt_C_over_fbar_[cell_counter-2]*implicit_rhs[Total_Cells_+cell_counter-2])/dx[i] ; 

		three_outlet_flux_[i][1] *= implicit_relaxation*S_x_[face_counter-1];

		three_outlet_flux_[i][2] = outlet_coef * (implicit_rhs[2*Total_Cells_+cell_counter-1]+implicit_rhs[2*Total_Cells_+cell_counter-2])/2.0 - 1.0/sqrt_Cf_x_[face_counter-1] * (sqrt_C_over_fbar_[cell_counter-1]*implicit_rhs[2*Total_Cells_+cell_counter-1] - sqrt_C_over_fbar_[cell_counter-2]*implicit_rhs[2*Total_Cells_+cell_counter-2])/dx[i];

		three_outlet_flux_[i][2] *= implicit_relaxation*S_x_[face_counter-1];

    }

    	solver.find_reservoir_dc(Channel_Inlet_Flux_,Channel_Outlet_Flux_, three_inlet_flux_, three_outlet_flux_,  Reservoir_Pressure_Type, Reservoir_Potential_Type, Connectivity, Connecting_Channel, Reservoir_Volume, dt, reservoir_dc);

    	cell_counter = 0;
    	for(int i=0; i<Num_Channel_; i++){
    	
			for(int j=0; j<Channel_Num_Cells_[i]+2; j++){
			
	    		dC_bar[cell_counter] = implicit_rhs[cell_counter] + 
	      		reservoir_dc[Channel_End_Reservoir[i][0]]*implicit_rhs[Total_Cells_+cell_counter] +
	      		reservoir_dc[Channel_End_Reservoir[i][1]]*implicit_rhs[2*Total_Cells_+cell_counter];
	      
	    		cell_counter++;
	  		}
      	}

    return;
}





