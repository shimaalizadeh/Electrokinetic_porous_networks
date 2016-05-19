#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<stdio.h>
#include <cmath>
using namespace std;


void Mesh_Generation(int Num_Channel, int Num_Reservoir, int *Ns, int *Nt, int *Nz, double *ds, double *dn, double *dz, double *theta, double *x0, double *y0, double *z0);

void write(int Num_Channel, int Num_Reservoir, int *Ns, int *Nt, int *Nz, double *C_bar, double *Reservoir_C_bar, double *C0, double *Reservoir_C0, double *Pressure, double *Reservoir_Pressure, double *Potential,double *Reservoir_Potential,  double time);


////////Function Implementations/////////
////////
void Mesh_Generation(int Num_Channel, int Num_Reservoir, int *Ns, int *Nt, int *Nz, 
double *ds, double *dn, double *dz, double *theta, double *x0, double *y0, double *z0){

  ////Creating Mesh file
  ofstream Myfile_Mesh;

  int domain_num=Num_Channel+Num_Reservoir;
  double X,Y,Z;  

  Myfile_Mesh.open("Mesh.bin", ios::out | ios::binary);
  Myfile_Mesh.write((char*) & domain_num, sizeof(int));

 for(int region=0; region<domain_num; region++)
   {
     //Note:&N[region] is equal to N+region where N is pointer array!
     Myfile_Mesh.write((char*) &Ns[region], sizeof(int));
     Myfile_Mesh.write((char*) &Nt[region], sizeof(int));
     Myfile_Mesh.write((char*) &Nz[region], sizeof(int));
   }

for(int region=0; region<domain_num; region++)
{
				
			for (int k=0;k<Nz[region];k++) {
				for (int j=0;j<Nt[region];j++) {
					for (int i=0;i<Ns[region];i++) {

						X=x0[region]+(i+0.5)*ds[region]*cos(theta[region]) +
							
							-(j+0.5)*dn[region]*sin(theta[region]);

						Myfile_Mesh.write((char*) & X, sizeof(double));
					}
				}
			}

			for (int k=0;k<Nz[region];k++) {
				for (int j=0;j<Nt[region];j++) {
					for (int i=0;i<Ns[region];i++) {

						Y=y0[region]+(j+1./2.)*dn[region]*cos(theta[region]) + 
							
							(i+1./2.)*ds[region]*sin(theta[region]);

						Myfile_Mesh.write((char*) & Y, sizeof(double));
					}
				}   
			}           
			
			for (int k=0;k<Nz[region];k++) {
				for (int j=0;j<Nt[region];j++) {
					for (int i=0;i<Ns[region];i++) {

					        Z=z0[region]+(k+1./2.)*dz[region];

						Myfile_Mesh.write((char*) & Z, sizeof(double));
					}
				}   
			} 

}

Myfile_Mesh.close();

return;
}

void write(int Num_Channel, int Num_Reservoir, int *Ns, int *Nt, int *Nz, double *C_bar, 
double *Reservoir_C_bar, double *C0, double *Reservoir_C0, double *Pressure, 
double *Reservoir_Pressure, double *Potential,double *Reservoir_Potential,  int time)
{
  int domain_num=Num_Channel+Num_Reservoir;

  int variable_num=4;

  double C;
  double P, Phi;

  ////writing filename
  ofstream Myfile_Data_Tec;
  string filename_base_Tec="Tecplot";
  ostringstream filename_out_Data_Tec;

  if(time < 10)
        filename_out_Data_Tec<<filename_base_Tec<<"_00000"<<time<<".bin";
  else if(time < 100)
        filename_out_Data_Tec<<filename_base_Tec<<"_0000"<<time<<".bin";
  else if(time <1000)
        filename_out_Data_Tec<<filename_base_Tec<<"_000"<<time<<".bin";
  else if(time < 10000)
        filename_out_Data_Tec<<filename_base_Tec<<"_00"<<time<<".bin";
  else if(time < 100000)
        filename_out_Data_Tec<<filename_base_Tec<<"_0"<<time<<".bin";
  else
  	filename_out_Data_Tec<<filename_base_Tec<<"_"<<time<<".bin";

  string filename_Tec=filename_out_Data_Tec.str();
	

 Myfile_Data_Tec.open(filename_Tec.c_str(), ios::out|ios::binary);
 Myfile_Data_Tec.write((char*) & domain_num, sizeof(int));			

  for(int region=0; region<domain_num; region++)
  {
	  Myfile_Data_Tec.write((char*) & Ns[region], sizeof(int));
	  Myfile_Data_Tec.write((char*) & Nt[region], sizeof(int));
	  Myfile_Data_Tec.write((char*) & Nz[region], sizeof(int));
	  Myfile_Data_Tec.write((char*) & variable_num, sizeof(int));
  }

  int counter=1;

  for(int region=0; region<domain_num; region++)
    {
      if(region< Num_Channel)
	{	


	  //C_bar		
	  for (int k=0;k<Nz[region];k++) {
		   for (int j=0;j<Nt[region];j++) {
		     for (int i=0;i<Ns[region];i++) {

		       C=C_bar[counter+i];

		       Myfile_Data_Tec.write((char*) & C, sizeof(double));
		     }
		   }
		 }
	  //C0
	  for (int k=0;k<Nz[region];k++) {
	    for (int j=0;j<Nt[region];j++) {
	      for (int i=0;i<Ns[region];i++) {

		C= C0[counter+i];

		Myfile_Data_Tec.write((char*) & C, sizeof(double));
	      }
	    }
	  }
	  
	  //Pressure
	  for (int k=0;k<Nz[region];k++) {
	    for (int j=0;j<Nt[region];j++) {
	      for (int i=0;i<Ns[region];i++) {

		P=Pressure[counter+i];

		Myfile_Data_Tec.write((char*) & P, sizeof(double));
	      }
	    }
	  }

	  //Potential
	  for (int k=0;k<Nz[region];k++) {
	    for (int j=0;j<Nt[region];j++) {
	      for (int i=0;i<Ns[region];i++) {

		Phi=Potential[counter+i];

		Myfile_Data_Tec.write((char*) & Phi, sizeof(double));
	      }
	    }
	  }
	  
	  counter+=Ns[region]+2;
	}//end of if

      else{

	//C_bar			
	for (int k=0;k<Nz[region];k++) {
	  for (int j=0;j<Nt[region];j++) {
	    for (int i=0;i<Ns[region];i++) {

	      C=Reservoir_C_bar[region-Num_Channel];

	      Myfile_Data_Tec.write((char*) & C, sizeof(double));
	    }
	  }
	}

	//C0
	for (int k=0;k<Nz[region];k++) {
	  for (int j=0;j<Nt[region];j++) {
	    for (int i=0;i<Ns[region];i++) {

	      C=Reservoir_C0[region-Num_Channel];

	      Myfile_Data_Tec.write((char*) & C, sizeof(double));
	    }
	  }
	}
	
	//Pressure
	for (int k=0;k<Nz[region];k++) {
	  for (int j=0;j<Nt[region];j++) {
	    for (int i=0;i<Ns[region];i++) {

	      P=Reservoir_Pressure[region-Num_Channel];

	      Myfile_Data_Tec.write((char*) & P, sizeof(double));
	    }
	  }
	}

	//Potential
	for (int k=0;k<Nz[region];k++) {
	  for (int j=0;j<Nt[region];j++) {
	    for (int i=0;i<Ns[region];i++) {

	      Phi=Reservoir_Potential[region-Num_Channel];

	      Myfile_Data_Tec.write((char*) & Phi, sizeof(double));
	    }
	  }
	}
	
      }//end of else
    }//end of for loop
			
  Myfile_Data_Tec.close();

return;
}
//////////
//////////End of function Implementation///////////


