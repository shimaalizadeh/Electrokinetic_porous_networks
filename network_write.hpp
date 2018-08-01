#ifndef NETWORK_WRITE_HPP
#define NETWORK_WRITE_HPP

#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<stdio.h>
#include <cmath>
using namespace std;


void Mesh_Generation(int Num_Channel, int Num_Reservoir, int *Ns, int *Nt, int *Nz, double *ds, double *dn, double *dz, double *theta, double *x0, double *y0, double *z0);

void write(int Num_Channel, int Num_Reservoir, int *Ns, int *Nt, int *Nz, double *C_bar, double *Reservoir_C_bar, double *C0, double *Reservoir_C0, double *Pressure, double *Reservoir_Pressure, double *Potential,double *Reservoir_Potential, double* U, int num_h_pores, int time);

void write_current_flow(int Num_channel, double *Q, double *I, double *S);


/*********IMPLEMENTATIONS*********/

void Mesh_Generation(int Num_Channel, int Num_Reservoir, int *Ns, int *Nt, int *Nz, 
double *ds, double *dn, double *dz, double *theta, double *x0, double *y0, double *z0){

  ////Creating Mesh file
  ofstream Myfile_Mesh;

  int domain_num=Num_Channel;
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
double *Reservoir_Pressure, double *Potential,double *Reservoir_Potential, double *U, int num_h_pores, int time)
{
  int domain_num = Num_Channel;

  int variable_num = 4 + 3;

  double C;
  double P, Phi;
  double u, v, w;

  ////writing filename
  ofstream Myfile_Data_Tec;
  string filename_base_Tec="Tecplot";
  ostringstream filename_out_Data_Tec;

  if(time < 10)
        filename_out_Data_Tec<<filename_base_Tec<<"_000"<<time<<".bin";
  else if(time < 100)
        filename_out_Data_Tec<<filename_base_Tec<<"_00"<<time<<".bin";
  else if(time <1000)
        filename_out_Data_Tec<<filename_base_Tec<<"_0"<<time<<".bin";
  else
        filename_out_Data_Tec<<filename_base_Tec<<"_"<<time<<".bin";

  string filename_Tec = filename_out_Data_Tec.str();

 Myfile_Data_Tec.open(filename_Tec.c_str(), ios::out|ios::binary);
 Myfile_Data_Tec.write((char*) & domain_num, sizeof(int));			

  for(int region=0; region<domain_num; region++)
  {
	  Myfile_Data_Tec.write((char*) & Ns[region], sizeof(int));
	  Myfile_Data_Tec.write((char*) & Nt[region], sizeof(int));
	  Myfile_Data_Tec.write((char*) & Nz[region], sizeof(int));
	  Myfile_Data_Tec.write((char*) & variable_num, sizeof(int));
  }

  int counter = 1;
  int face_counter = 0;

  for(int region=0; region<domain_num; region++)
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

	  //U
	  for (int k=0;k<Nz[region];k++) {
            for (int j=0;j<Nt[region];j++) {
              for (int i=0;i<Ns[region];i++) {

		if(region < num_h_pores)
			u = U[face_counter+i+1];
		else
			u = 0.0;

                Myfile_Data_Tec.write((char*) &u, sizeof(double));
              }
            }
          }
	  
	  //V
	  for (int k=0;k<Nz[region];k++) {
            for (int j=0;j<Nt[region];j++) {
              for (int i=0;i<Ns[region];i++) {

		if(region < num_h_pores)
			v = 0.0;
		else
                	v = U[face_counter+i+1];

                Myfile_Data_Tec.write((char*) &v, sizeof(double));
              }
            }
          }

	  //W
	  for (int k=0;k<Nz[region];k++) {
            for (int j=0;j<Nt[region];j++) {
              for (int i=0;i<Ns[region];i++) {

                w = 0.0;

                Myfile_Data_Tec.write((char*) &w, sizeof(double));
              }
            }
          }

	  counter += Ns[region]+2;
	  face_counter += Ns[region]+1;
	}
 Myfile_Data_Tec.close();

return;
}


void write_current_flow(int Num_channel, double *Q, double *I, double *S){

  ofstream fwrite;
  int numchannel = Num_channel;

  fwrite.open("current_flow_area.bin", ios::out | ios::binary);
  fwrite.write((char*) &numchannel, sizeof(int));

  fwrite.write((char*) Q, numchannel*sizeof(double));
  fwrite.write((char*) I, numchannel*sizeof(double));
  fwrite.write((char*) S, numchannel*sizeof(double)); 

  return;
}

#endif
