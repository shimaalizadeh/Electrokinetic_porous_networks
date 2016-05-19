#include"network_info.hpp"
#include <stdlib.h>


void Network::setup(const char* filename_1, int Num_Channel_, int Num_Reservoir_, int restart_, double right_voltage){
    
  //set Pe and gp_bar
  set_electrochemical_parameters();
  circle_gp_bar = -1./2.;  //for cicular pore with R=2
  slit_gp_bar = -1./3.; //for slit pore

  Num_Channel = Num_Channel_;
 
  Num_Reservoir = Num_Reservoir_;
  
  restart = restart_;

  Initial_Memory_Allocation(Num_Channel, Num_Reservoir);
  
  Read_filename_1(filename_1, right_voltage);
  
  //I want to force zero value for P, Phi, k and dC_bar at the beginning
  set_to_zero(Pressure, Total_Cells);
  set_to_zero(Potential, Total_Cells);
  set_to_zero(dC_bar, Total_Cells);
  
}

void Network::reset(const char* filename_1, double right_voltage) {

	restart = 0;
	time = 0.;
	reset_from_file(filename_1, right_voltage);
	
	return;
}

////Destructor
Network::~Network(){

  ////free memory
  delete [] sigma_val;
    delete [] lambda_val;
    delete [] s_val;
    delete [] c0_val;
    
  delete [] Channel_type;
  delete [] Reservoir_type;
  delete [] Channel_Num_Cells;
  delete [] Reservoir_Pressure_Type;
  delete [] Reservoir_Potential_Type;
  delete [] Reservoir_Pressure_Value; 
  delete [] Reservoir_Potential_Value;
  delete [] Reservoir_C_bar;
  delete [] Reservoir_C0;
  delete [] Reservoir_Sigma;
  delete [] Reservoir_Lambda;
  delete [] Reservoir_Pressure;
  delete [] Reservoir_Potential;
  delete [] Reservoir_Volume;
  delete [] dx;
  delete [] C_bar;
  delete [] C0;
  delete [] Cs;
  delete [] sigma_star;
  delete [] lambda;
  delete [] lambda_for_table;
  delete [] nondimensional_Sp;
  delete [] Pressure;
  delete [] Potential;
  delete [] U;
  delete [] dPdx;
  delete [] dMudx;
  delete [] Art_Diff;
  delete [] Q1;
  delete [] I1;
  delete [] Q2;
  delete [] I2;
  delete [] Q3;
  delete [] I3;
  delete [] Q;
  delete [] I;
  delete [] dC_bar;
  delete [] Ns;
  delete [] Nt;
  delete [] Nz;
  delete [] ds;
  delete [] dn;
  delete [] dz;
  delete [] x0;
  delete [] y0;
  delete [] z0;
  delete [] theta;

  for(int i=0;i<Num_Channel;i++)
    {
      delete [] Channel_End_Reservoir[i];
    }
  delete [] Channel_End_Reservoir;

  for(int i=0;i<Num_Reservoir;i++)
    {
      delete [] Connectivity[i];
      delete [] Connecting_Channel[i];
    }
  delete [] Connectivity;
  delete [] Connecting_Channel;

}


void Network::Initial_Memory_Allocation(int  Num_Channel, int Num_Reservoir){

  sigma_val = new double [Num_Channel];
  lambda_val = new double [Num_Channel];
  s_val = new double [Num_Channel];
  c0_val = new double [Num_Channel];

  Channel_type = new int[Num_Channel];
 
  Reservoir_type = new int[Num_Reservoir];
  
  Ns=new int[Num_Channel+Num_Reservoir];

  Nt=new int[Num_Channel+Num_Reservoir];

  Nz=new int[Num_Channel+Num_Reservoir];

  ds=new double[Num_Channel+Num_Reservoir];
  
  dn=new double[Num_Channel+Num_Reservoir];
  
  dz=new double[Num_Channel+Num_Reservoir]; 
  
  x0=new double[Num_Channel+Num_Reservoir];
  
  y0=new double[Num_Channel+Num_Reservoir];
  
  z0=new double[Num_Channel+Num_Reservoir];
  
  theta=new double[Num_Channel+Num_Reservoir];


  Channel_End_Reservoir=new int*[Num_Channel];

  for(int i=0;i<Num_Channel;i++)
    {
      Channel_End_Reservoir[i]=new int[2];
    }

  Channel_Num_Cells=new int[Num_Channel];

  dx=new double[Num_Channel];

  Reservoir_Pressure_Type=new bool[Num_Reservoir];

  Reservoir_Potential_Type=new bool[Num_Reservoir];

  Reservoir_Pressure_Value=new double[Num_Reservoir];
  
  Reservoir_Potential_Value=new double[Num_Reservoir];

  Reservoir_C_bar=new double[Num_Reservoir];

  Reservoir_C0=new double[Num_Reservoir];

  Reservoir_Sigma=new double[Num_Reservoir];

  Reservoir_Lambda=new double[Num_Reservoir];

  Reservoir_Pressure=new double[Num_Reservoir];

  Reservoir_Potential=new double[Num_Reservoir];

  Reservoir_Volume=new double[Num_Reservoir];
  
  Connectivity= new int*[Num_Reservoir];
  
  Connecting_Channel=new int*[Num_Reservoir];

  for(int i=0;i<Num_Reservoir;i++)
    {
      Connectivity[i]=new int[Num_Reservoir];

      Connecting_Channel[i]=new int[Num_Reservoir];
    }

return;
}

void Network::Read_filename_1(const char* filename_1, double right_voltage){

  ifstream read_1(filename_1);
  double S_ref;
  

  // Reading if artificial diffsivity is needed
  read_1>>Art_Diff_Switch;
  
  if(restart == 0)
  {  
    // Reading Channel properties
	for(int i=0; i<Num_Channel; i++)
	{
	  read_1>>Channel_type[i]>>Channel_End_Reservoir[i][0]>>Channel_End_Reservoir[i][1]>>dx[i]>>Channel_Num_Cells[i]>>
	  sigma_val[i]>>lambda_val[i]>>s_val[i]>>c0_val[i];
	  
	  dx[i]=dx[i]/Channel_Num_Cells[i]; 
	 }

    // Reading Reservoir properties
   	for(int i=0;i<Num_Reservoir;i++)
	{
		  read_1>>Reservoir_type[i]>>Reservoir_Pressure_Type[i]>>Reservoir_Potential_Type[i]>>Reservoir_Pressure_Value[i]>>Reservoir_Potential_Value[i]>>Reservoir_C0[i]>>Reservoir_Sigma[i]>>Reservoir_Lambda[i]>>Reservoir_Volume[i];

			
		  Reservoir_Potential_Value[i] = Reservoir_Potential_Value[i]; //EDITTED
		  Reservoir_Pressure_Value[i] = -2*Reservoir_C0[i] + Reservoir_Pressure_Value[i];
		
		  //Reservoir_Potential_Value[i] = 0.;
		  //Reservoir_Pressure_Value[i] = 0.;	
	}
	
		
      Total_Cells = 0;
      Total_Faces = 0;

  	  for(int i=0;i<Num_Channel;i++)
      {
      		Total_Cells += Channel_Num_Cells[i];
      		Total_Faces += Channel_Num_Cells[i];
      }
  		//for each channel we have two reservoirs at the ends:
  		Total_Cells += Num_Channel*2;
  		Total_Faces += Num_Channel;
  
   		Memory_Allocation();
   		
   		 cell_counter=0;
   		 for(int i=0; i<Num_Channel; i++){
   		
   			for(int j=1; j<Channel_Num_Cells[i]+1; j++){
   			
	  			sigma_star[cell_counter+j] = sigma_val[i];
	  		
	  			lambda[cell_counter+j] = lambda_val[i];
	  		
	  			nondimensional_Sp[cell_counter+j] = s_val[i];
	  		
	  			C0[cell_counter+j] = c0_val[i];

	  			Cs[cell_counter+j] = 2*lambda[cell_counter+j]*lambda[cell_counter+j]*(abs(sigma_star[cell_counter+j]));
	  		
          		lambda_for_table[cell_counter+j] = lambda[cell_counter+j];

	  			// initial value of Mu
	  			Potential[cell_counter+j] = log(C0[cell_counter+j]) + Potential[cell_counter+j];

				// initial P_tot
				Pressure[cell_counter+j] = -2*C0[cell_counter+j] + Pressure[cell_counter+j];
         
	   		}//End of for(int j...)
	   		
	   		// setting end reservoirs
    		//Note: lambda and S are set based on adjacent cell and they may not
    		// include the correct values. The correct value of reservoir lambda 
    		//is in lambda_for_table vector.
       
    		sigma_star[cell_counter] = Reservoir_Sigma[Channel_End_Reservoir[i][0]];

    		lambda_for_table[cell_counter] = Reservoir_Lambda[Channel_End_Reservoir[i][0]];

    		lambda[cell_counter] = lambda[cell_counter+1];

      		nondimensional_Sp[cell_counter] = nondimensional_Sp[cell_counter+1];

      		Cs[cell_counter] = 2*Reservoir_Lambda[Channel_End_Reservoir[i][0]]*Reservoir_Lambda[Channel_End_Reservoir[i][0]]*abs(Reservoir_Sigma[Channel_End_Reservoir[i][0]]);

      		C0[cell_counter] = Reservoir_C0[Channel_End_Reservoir[i][0]];

      		Potential[cell_counter] = Reservoir_Potential_Value[Channel_End_Reservoir[i][0]];

      		Pressure[cell_counter] = Reservoir_Pressure_Value[Channel_End_Reservoir[i][0]];

      		cell_counter += Channel_Num_Cells[i]+2;

      		sigma_star[cell_counter-1] = Reservoir_Sigma[Channel_End_Reservoir[i][1]];

      		lambda_for_table[cell_counter-1] = Reservoir_Lambda[Channel_End_Reservoir[i][1]];

      		lambda[cell_counter-1] = lambda[cell_counter-2];

      		nondimensional_Sp[cell_counter-1] = nondimensional_Sp[cell_counter-2];

      		Cs[cell_counter-1] = 2*Reservoir_Lambda[Channel_End_Reservoir[i][1]]*Reservoir_Lambda[Channel_End_Reservoir[i][1]]*abs(Reservoir_Sigma[Channel_End_Reservoir[i][1]]);

      		C0[cell_counter-1] = Reservoir_C0[Channel_End_Reservoir[i][1]];

      		Potential[cell_counter-1] = Reservoir_Potential_Value[Channel_End_Reservoir[i][1]];

      		Pressure[cell_counter-1] = Reservoir_Pressure_Value[Channel_End_Reservoir[i][1]];
      	}
	
    }
    else
      {
		// Reading Channel properties
		for(int i=0; i<Num_Channel; i++)
		{
	  		read_1>>Channel_type[i]>>Channel_End_Reservoir[i][0]>>Channel_End_Reservoir[i][1]>>dx[i]>>Channel_Num_Cells[i]>>
	  		sigma_val[i]>>lambda_val[i]>>s_val[i]>>c0_val[i];	
	  		
        }

		// Reading Reservoir properties
		for(int i=0; i<Num_Reservoir; i++)
	  	{
	    	read_1>>Reservoir_type[i]>>Reservoir_Pressure_Type[i]>>Reservoir_Potential_Type[i]>>Reservoir_Pressure_Value[i]>>Reservoir_Potential_Value[i]>>Reservoir_C_bar[i]>>Reservoir_Sigma[i]>>Reservoir_Lambda[i]>>Reservoir_Volume[i];

	  	}
	  	
	Total_Cells = 0;
  	Total_Faces = 0;

  	for(int i=0; i<Num_Channel; i++)
    {
      		Total_Cells += Channel_Num_Cells[i];
      		Total_Faces += Channel_Num_Cells[i];
    }
  	
  	//for each channel we have two reservoirs at the ends:
  	Total_Cells += Num_Channel*2;
  	Total_Faces += Num_Channel;
  
   	Memory_Allocation();
	  	
	cell_counter = 0;
    
    for(int i=0; i<Num_Channel; i++){
    
		for(int j=1; j<Channel_Num_Cells[i]+1; j++){
		
	  			sigma_star[cell_counter+j] = sigma_val[i];
	  		
	  			lambda[cell_counter+j] = lambda_val[i];
	  		
	  			nondimensional_Sp[cell_counter+j] = s_val[i];
	  		
	  			C_bar[cell_counter+j] = c0_val[i];

	  			Cs[cell_counter+j] = 2*lambda[cell_counter+j]*lambda[cell_counter+j]*(abs(sigma_star[cell_counter+j]));
	  		
          		lambda_for_table[cell_counter+j] = lambda[cell_counter+j];

	   	}//End of for(int j...)	
	   	
    	sigma_star[cell_counter] = Reservoir_Sigma[Channel_End_Reservoir[i][0]];

      	lambda_for_table[cell_counter] = Reservoir_Lambda[Channel_End_Reservoir[i][0]];

      	lambda[cell_counter] = lambda[cell_counter+1];

      	nondimensional_Sp[cell_counter] = nondimensional_Sp[cell_counter+1];

      	Cs[cell_counter] = 2*Reservoir_Lambda[Channel_End_Reservoir[i][0]]*Reservoir_Lambda[Channel_End_Reservoir[i][0]]*abs(Reservoir_Sigma[Channel_End_Reservoir[i][0]]);

      	C_bar[cell_counter] = Reservoir_C_bar[Channel_End_Reservoir[i][0]];

      	Potential[cell_counter] = Reservoir_Potential_Value[Channel_End_Reservoir[i][0]];

      	Pressure[cell_counter] = Reservoir_Pressure_Value[Channel_End_Reservoir[i][0]];

      	cell_counter += Channel_Num_Cells[i]+2;

      	sigma_star[cell_counter-1] = Reservoir_Sigma[Channel_End_Reservoir[i][1]];

      	lambda_for_table[cell_counter-1] = Reservoir_Lambda[Channel_End_Reservoir[i][1]];

      	lambda[cell_counter-1] = lambda[cell_counter-2];

      	nondimensional_Sp[cell_counter-1] = nondimensional_Sp[cell_counter-2];

      	Cs[cell_counter-1] = 2*Reservoir_Lambda[Channel_End_Reservoir[i][1]]*Reservoir_Lambda[Channel_End_Reservoir[i][1]]*abs(Reservoir_Sigma[Channel_End_Reservoir[i][1]]);

      	C_bar[cell_counter-1] = Reservoir_C_bar[Channel_End_Reservoir[i][1]];

      	Potential[cell_counter-1] = Reservoir_Potential_Value[Channel_End_Reservoir[i][1]];

      	Pressure[cell_counter-1] = Reservoir_Pressure_Value[Channel_End_Reservoir[i][1]];
    }
}
	
    // Reading connectivity matrix:
	int row, column;
	
	for(int i=0; i<Num_Reservoir; i++)
	{
		for(int j=0; j<Num_Reservoir; j++)
		{
			Connectivity[i][j] = 0;

			Connecting_Channel[i][j] = 0;
		}
	}

	for(int i=0; i<Num_Channel; i++)
	{	
	  read_1>>row>>column;

	  Connectivity[row][column] = 1;

	  Connectivity[column][row] = -1;

	  read_1>>Connecting_Channel[row][column];
			
	  Connecting_Channel[column][row] = Connecting_Channel[row][column];

	}
    
    //geometry parameters
	if(restart == 0)
	  {
	    for(int i=0;i<Num_Channel+Num_Reservoir;i++)
	      {
			read_1>>Ns[i]>>Nt[i]>>Nz[i]>>ds[i]>>dn[i]>>dz[i]>>x0[i]>>y0[i]>>z0[i]>>theta[i];
      
			ds[i] = ds[i]/(Ns[i]);

			dn[i] = dn[i]/(Nt[i]);

			dz[i] = dz[i]/(Nz[i]);
	      }
          }
	else
	  {
	    for(int i=0; i<Num_Channel+Num_Reservoir; i++)
	      {
			read_1>>Ns[i]>>Nt[i]>>Nz[i]>>ds[i]>>dn[i]>>dz[i]>>x0[i]>>y0[i]>>z0[i]>>theta[i];
	      }
	  }
	
	S_ref = nondimensional_Sp[1];
  	for(int i=0; i<Total_Cells; i++)
    {
      nondimensional_Sp[i]=nondimensional_Sp[i]/S_ref;
    }
    
	read_1.close();

  return;
}

void Network::reset_from_file(const char* filename_1, double right_voltage){

  ifstream read_1(filename_1);
  double S_ref;
  

  // Reading if artificial diffsivity is needed
  read_1>>Art_Diff_Switch;
  
  if(restart == 0)
  {  
    // Reading Channel properties
	for(int i=0; i<Num_Channel; i++)
	{
	  read_1>>Channel_type[i]>>Channel_End_Reservoir[i][0]>>Channel_End_Reservoir[i][1]>>dx[i]>>Channel_Num_Cells[i]>>
	  sigma_val[i]>>lambda_val[i]>>s_val[i]>>c0_val[i];
	  
	  dx[i]=dx[i]/Channel_Num_Cells[i]; 
	 }

    // Reading Reservoir properties
   	for(int i=0;i<Num_Reservoir;i++)
	{
		  read_1>>Reservoir_type[i]>>Reservoir_Pressure_Type[i]>>Reservoir_Potential_Type[i]>>Reservoir_Pressure_Value[i]>>Reservoir_Potential_Value[i]>>Reservoir_C0[i]>>Reservoir_Sigma[i]>>Reservoir_Lambda[i]>>Reservoir_Volume[i];

		  Reservoir_Potential_Value[i] = log(Reservoir_C0[i]) + Reservoir_Potential_Value[i];
		  Reservoir_Pressure_Value[i] = -2*Reservoir_C0[i] + Reservoir_Pressure_Value[i];
			
	}
	
	//update right reservoir potential///////////////////////////////////
	//Reservoir_Potential_Value[2] = right_voltage + log(Reservoir_C0[2]);
	/////////////////////////////////////////////////////////////////////
		
   		 cell_counter=0;
   		 for(int i=0; i<Num_Channel; i++){
   		
   			for(int j=1; j<Channel_Num_Cells[i]+1; j++){
   			
	  			sigma_star[cell_counter+j] = sigma_val[i];
	  		
	  			lambda[cell_counter+j] = lambda_val[i];
	  		
	  			nondimensional_Sp[cell_counter+j] = s_val[i];
	  		
	  			C0[cell_counter+j] = c0_val[i];

	  			Cs[cell_counter+j] = 2*lambda[cell_counter+j]*lambda[cell_counter+j]*(abs(sigma_star[cell_counter+j]));
	  		
          		lambda_for_table[cell_counter+j] = lambda[cell_counter+j];

	  			// initial value of Mu
	  			Potential[cell_counter+j] = log(C0[cell_counter+j]) + Potential[cell_counter+j];

				// initial P_tot
				Pressure[cell_counter+j] = -2*C0[cell_counter+j] + Pressure[cell_counter+j];
         
	   		}//End of for(int j...)
	   		
	   		// setting end reservoirs
    		//Note: lambda and S are set based on adjacent cell and they may not
    		// include the correct values. The correct value of reservoir lambda 
    		//is in lambda_for_table vector.
       
    		sigma_star[cell_counter] = Reservoir_Sigma[Channel_End_Reservoir[i][0]];

    		lambda_for_table[cell_counter] = Reservoir_Lambda[Channel_End_Reservoir[i][0]];

    		lambda[cell_counter] = lambda[cell_counter+1];

      		nondimensional_Sp[cell_counter] = nondimensional_Sp[cell_counter+1];

      		Cs[cell_counter] = 2*Reservoir_Lambda[Channel_End_Reservoir[i][0]]*Reservoir_Lambda[Channel_End_Reservoir[i][0]]*abs(Reservoir_Sigma[Channel_End_Reservoir[i][0]]);

      		C0[cell_counter] = Reservoir_C0[Channel_End_Reservoir[i][0]];

      		Potential[cell_counter] = Reservoir_Potential_Value[Channel_End_Reservoir[i][0]];

      		Pressure[cell_counter] = Reservoir_Pressure_Value[Channel_End_Reservoir[i][0]];

      		cell_counter += Channel_Num_Cells[i]+2;

      		sigma_star[cell_counter-1] = Reservoir_Sigma[Channel_End_Reservoir[i][1]];

      		lambda_for_table[cell_counter-1] = Reservoir_Lambda[Channel_End_Reservoir[i][1]];

      		lambda[cell_counter-1] = lambda[cell_counter-2];

      		nondimensional_Sp[cell_counter-1] = nondimensional_Sp[cell_counter-2];

      		Cs[cell_counter-1] = 2*Reservoir_Lambda[Channel_End_Reservoir[i][1]]*Reservoir_Lambda[Channel_End_Reservoir[i][1]]*abs(Reservoir_Sigma[Channel_End_Reservoir[i][1]]);

      		C0[cell_counter-1] = Reservoir_C0[Channel_End_Reservoir[i][1]];

      		Potential[cell_counter-1] = Reservoir_Potential_Value[Channel_End_Reservoir[i][1]];

      		Pressure[cell_counter-1] = Reservoir_Pressure_Value[Channel_End_Reservoir[i][1]];
      	}
	
    }
    else
      {
		// Reading Channel properties
		for(int i=0; i<Num_Channel; i++)
		{
	  		read_1>>Channel_type[i]>>Channel_End_Reservoir[i][0]>>Channel_End_Reservoir[i][1]>>dx[i]>>Channel_Num_Cells[i]>>
	  		sigma_val[i]>>lambda_val[i]>>s_val[i]>>c0_val[i];	
	  		
        }

		// Reading Reservoir properties
		for(int i=0; i<Num_Reservoir; i++)
	  	{
	    	read_1>>Reservoir_type[i]>>Reservoir_Pressure_Type[i]>>Reservoir_Potential_Type[i]>>Reservoir_Pressure_Value[i]>>Reservoir_Potential_Value[i]>>Reservoir_C_bar[i]>>Reservoir_Sigma[i]>>Reservoir_Lambda[i]>>Reservoir_Volume[i];

	  	}
	  	
	cell_counter = 0;
    
    for(int i=0; i<Num_Channel; i++){
    
		for(int j=1; j<Channel_Num_Cells[i]+1; j++){
		
	  			sigma_star[cell_counter+j] = sigma_val[i];
	  		
	  			lambda[cell_counter+j] = lambda_val[i];
	  		
	  			nondimensional_Sp[cell_counter+j] = s_val[i];
	  		
	  			C_bar[cell_counter+j] = c0_val[i];

	  			Cs[cell_counter+j] = 2*lambda[cell_counter+j]*lambda[cell_counter+j]*(abs(sigma_star[cell_counter+j]));
	  		
          		lambda_for_table[cell_counter+j] = lambda[cell_counter+j];

	   	}//End of for(int j...)	
	   	
    	sigma_star[cell_counter] = Reservoir_Sigma[Channel_End_Reservoir[i][0]];

      	lambda_for_table[cell_counter] = Reservoir_Lambda[Channel_End_Reservoir[i][0]];

      	lambda[cell_counter] = lambda[cell_counter+1];

      	nondimensional_Sp[cell_counter] = nondimensional_Sp[cell_counter+1];

      	Cs[cell_counter] = 2*Reservoir_Lambda[Channel_End_Reservoir[i][0]]*Reservoir_Lambda[Channel_End_Reservoir[i][0]]*abs(Reservoir_Sigma[Channel_End_Reservoir[i][0]]);

      	C_bar[cell_counter] = Reservoir_C_bar[Channel_End_Reservoir[i][0]];

      	Potential[cell_counter] = Reservoir_Potential_Value[Channel_End_Reservoir[i][0]];

      	Pressure[cell_counter] = Reservoir_Pressure_Value[Channel_End_Reservoir[i][0]];

      	cell_counter += Channel_Num_Cells[i]+2;

      	sigma_star[cell_counter-1] = Reservoir_Sigma[Channel_End_Reservoir[i][1]];

      	lambda_for_table[cell_counter-1] = Reservoir_Lambda[Channel_End_Reservoir[i][1]];

      	lambda[cell_counter-1] = lambda[cell_counter-2];

      	nondimensional_Sp[cell_counter-1] = nondimensional_Sp[cell_counter-2];

      	Cs[cell_counter-1] = 2*Reservoir_Lambda[Channel_End_Reservoir[i][1]]*Reservoir_Lambda[Channel_End_Reservoir[i][1]]*abs(Reservoir_Sigma[Channel_End_Reservoir[i][1]]);

      	C_bar[cell_counter-1] = Reservoir_C_bar[Channel_End_Reservoir[i][1]];

      	Potential[cell_counter-1] = Reservoir_Potential_Value[Channel_End_Reservoir[i][1]];

      	Pressure[cell_counter-1] = Reservoir_Pressure_Value[Channel_End_Reservoir[i][1]];
    }
}
	
    // Reading connectivity matrix:
	int row, column;
	
	for(int i=0; i<Num_Reservoir; i++)
	{
		for(int j=0; j<Num_Reservoir; j++)
		{
			Connectivity[i][j] = 0;

			Connecting_Channel[i][j] = 0;
		}
	}

	for(int i=0; i<Num_Channel; i++)
	{	
	  read_1>>row>>column;

	  Connectivity[row][column] = 1;

	  Connectivity[column][row] = -1;

	  read_1>>Connecting_Channel[row][column];
			
	  Connecting_Channel[column][row] = Connecting_Channel[row][column];

	}
    
    //geometry parameters
	if(restart == 0)
	  {
	    for(int i=0;i<Num_Channel+Num_Reservoir;i++)
	      {
			read_1>>Ns[i]>>Nt[i]>>Nz[i]>>ds[i]>>dn[i]>>dz[i]>>x0[i]>>y0[i]>>z0[i]>>theta[i];
      
			ds[i] = ds[i]/(Ns[i]);

			dn[i] = dn[i]/(Nt[i]);

			dz[i] = dz[i]/(Nz[i]);
	      }
          }
	else
	  {
	    for(int i=0; i<Num_Channel+Num_Reservoir; i++)
	      {
			read_1>>Ns[i]>>Nt[i]>>Nz[i]>>ds[i]>>dn[i]>>dz[i]>>x0[i]>>y0[i]>>z0[i]>>theta[i];
	      }
	  }
	
	S_ref = nondimensional_Sp[1];
  	for(int i=0; i<Total_Cells; i++)
    {
      nondimensional_Sp[i]=nondimensional_Sp[i]/S_ref;
    }
    
	read_1.close();

  return;
}


void Network::Memory_Allocation(){

  sigma_star=new double[Total_Cells];

  lambda=new double[Total_Cells];

  lambda_for_table=new double[Total_Cells];

  nondimensional_Sp=new double[Total_Cells];

  C_bar=new double[Total_Cells];

  C0=new double[Total_Cells];

  Cs=new double[Total_Cells];

  Pressure=new double[Total_Cells];

  Potential=new double[Total_Cells];

  Q1=new double[Num_Channel];

  I1=new double[Num_Channel];

  Q2=new double[Num_Channel];

  I2=new double[Num_Channel];

  Q3=new double[Num_Channel];

  I3=new double[Num_Channel];

  Q=new double[Num_Channel];

  I=new double[Num_Channel];

  dC_bar=new double[Total_Cells];
    
  U=new double[Total_Faces];
  
  dPdx = new double[Total_Faces];
  
  dMudx = new double[Total_Faces];

  Art_Diff=new double[Total_Faces];

  return;
}

void Network::set_to_zero(double *vector, int size){
    for(int i=0; i<size; i++)
    {
        vector[i]=0.;
    }
    return;
}

void Network::set_electrochemical_parameters(){
    
    z=1;                 //valence
    e=1.602e-19;         //elementary charge(Coulombs!)
    Kb=1.38e-23;         //Boltzmann Constant
    epsilon0=8.85e-12;   //permitivity constant
    epsilon=80;          //dielectric relative permitivity
    mu=8.9e-4;           //dynamic viscosity(Pa.s)
    T=300;               //temperature(K)
    Dif=2e-9;            //mass Diffusivity (m^2/s),
    C_ref=0.001;         //reference concentration(Molar)
    NA=6.022e23; 
     
    lambda_ref=sqrt(epsilon*epsilon0*Kb*T/(2*C_ref*NA*1000*z*z*e*e));   
    Pe=epsilon0*epsilon/(mu*Dif)*(Kb*T/(z*e))*(Kb*T/(z*e));

    return;    
}

void Network::set_dt(){

  double *dt_=new double [Num_Channel];

  for(int i=0; i<Num_Channel; i++)   dt_[i]=dx[i]*dx[i]/4.;

  dt=dt_[0];

  for(int i=0; i<Num_Channel; i++)
    {
      if(dt>dt_[i]) dt=dt_[i];
    }

  delete [] dt_;
  return;
}

void Network::set_simulation_time(double time, double dt, double T_max, int period) {
	this->time = time;
	this->dt = dt;
	this->T_max = T_max;
	this->period = period;
	return;
}
void Network::store_reservoir_data(double *data, double*Reservoir_data){

  cell_counter=0;

  for(int i=0; i<Num_Channel; i++)
    {
      Reservoir_data[Channel_End_Reservoir[i][0]]=data[cell_counter];
      cell_counter+=Channel_Num_Cells[i]+2;
      Reservoir_data[Channel_End_Reservoir[i][1]]=data[cell_counter-1];
    }

 return;
}

void Network::Update_C_bar(int myid, double tt){

  for(int i=0; i<Total_Cells; i++)
    {
      C_bar[i]=C_bar[i]+dC_bar[i];
      dC_bar[i]=0.;

      if(abs(C_bar[i]) > 1000)
	{
	  cout<<"processor "<<myid<<" NAN concentration! at"<<i<<"  time="<<tt<<endl;
	  exit(EXIT_FAILURE);
	  }
    }
  
  ////Storing Reservoir C_bar
  store_reservoir_data(C_bar,Reservoir_C_bar); 

  return;
}

double * Network::Shock_Treatment(){

  double Diff_coeff;
  cell_counter=0;
  face_counter=0;

 

  for(int i=0; i<Num_Channel; i++)
    {
      for(int j=0; j<Channel_Num_Cells[i]+1; j++)
	{
	  Diff_coeff = 0.1*0.1/( (C_bar[cell_counter+j]/Cs[cell_counter+j]+C_bar[cell_counter+j+1]/Cs[cell_counter+j+1])/2.0 * (C_bar[cell_counter+j]/Cs[cell_counter+j]+C_bar[cell_counter+j+1]/Cs[cell_counter+j+1])/2.0 );

  Diff_coeff += 0.2*(C_bar[cell_counter+j+1]-C_bar[cell_counter+j])/( (C_bar[cell_counter+j+1]+C_bar[cell_counter+j])/2.0 )*(C_bar[cell_counter+j+1]-C_bar[cell_counter+j])/( (C_bar[cell_counter+j+1]+C_bar[cell_counter+j])/2.0 );

  Art_Diff[face_counter+j] = Diff_coeff/(1+Diff_coeff) * (dx[i]*I[i] +1)* heaviside(dx[i]*I[i] +1)* Art_Diff_Switch; 
  
  if(Art_Diff[face_counter+j] >= 0.0 )
	Art_Diff[face_counter+j] += 1.0;
  else
	Art_Diff[face_counter+j] = 1.0;

  //heaviside(-dx[i]*I[i]*sign((C_star[cell_counter+j+1]-C_star[cell_counter+j])/dx[j])-1)*Switch;
	}

      cell_counter+=Channel_Num_Cells[i]+2;
      face_counter+=Channel_Num_Cells[i]+1;

    }

  return(Art_Diff);
}

void Network::Write_filename_1(const char* filename_1)
{
  cout<<"writing network properties..."<<endl;
  cell_counter = 0;
  ofstream write_1(filename_1);
  write_1<<Art_Diff_Switch<<endl;

  for(int i=0; i<Num_Channel; i++)
    {
      write_1<<Channel_type[i]<<setw(5)<<Channel_End_Reservoir[i][0]<<setw(5)<<Channel_End_Reservoir[i][1]<<setw(10)<<dx[i]<<setw(10)<<Channel_Num_Cells[i]<<
      sigma_star[cell_counter+1]<<setw(10)<<lambda[cell_counter+1]<<setw(10)<<nondimensional_Sp[cell_counter+1]<<setw(10)<<
      C_bar[cell_counter+1]<<endl;
      
      cell_counter += Channel_Num_Cells[i]+2;
      
    }

  for(int i=0;i<Num_Reservoir;i++)
    {
      write_1<<Reservoir_type[i]<<setw(5)<<Reservoir_Pressure_Type[i]<<setw(5)<<Reservoir_Potential_Type[i]<<setw(15)<<Reservoir_Pressure[i]<<setw(20)<<Reservoir_Potential[i]<<setw(20)<<Reservoir_C_bar[i]<<setw(20)<<Reservoir_Sigma[i]<<setw(20)<<Reservoir_Lambda[i]<<setw(20)<<Reservoir_Volume[i]<<endl;                                                                                            
    }

  for(int i=0; i<Num_Channel; i++)
    {
      write_1<<Channel_End_Reservoir[i][0]<<setw(5)<<Channel_End_Reservoir[i][1]<<setw(5)<<i<<endl;
    }

  write_1<<endl;
  //geometry parameters
  for(int i=0;i<Num_Channel+Num_Reservoir;i++)
    {
      write_1<<Ns[i]<<setw(5)<<Nt[i]<<setw(5)<<Nz[i]<<setw(10)<<ds[i]<<setw(10)<<dn[i]<<setw(10)<<dz[i]<<setw(10)<<
	x0[i]<<setw(10)<<y0[i]<<setw(10)<<z0[i]<<setw(10)<<theta[i]<<endl;;
    }
     
  return;
}

double  Network::check_C_conservativity(void)
{
  // sum channel elements first
  double Total_C = 0.0;
  int cell_counter = 0;
  for(int i=0; i<Num_Channel; i++)
    {
      for(int j=1; j<Channel_Num_Cells[i]+1; j++)
	{
	  Total_C += C_bar[cell_counter+j]*dx[i]*nondimensional_Sp[cell_counter+j];
	}
      cell_counter += Channel_Num_Cells[i]+2;
    }

  // sum the reservoir concentrations
  for(int i=0; i<Num_Reservoir; i++)
    {
      Total_C += Reservoir_C_bar[i]*dx[0]*nondimensional_Sp[0];
    }

  return Total_C;
}

void Network::compute_gradients(void)
{
    cell_counter=0;
    face_counter=0;

    ////on faces
    for(int i=0; i<Num_Channel; i++)
    {
        for(int j=0; j<Channel_Num_Cells[i]+1; j++)
        {
		dPdx[face_counter] = (Pressure[cell_counter+j+1]-Pressure[cell_counter+j])/dx[i];
		
		dMudx[face_counter] = (Potential[cell_counter+j+1]-Potential[cell_counter+j])/dx[i];
		face_counter++;
	}

	cell_counter+=Channel_Num_Cells[i]+2;
    }
  return;
}


  
 
