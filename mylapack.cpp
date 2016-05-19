#include "mylapack.hpp"

void Lapack::setup(unsigned Num_Channel, unsigned Num_Reservoir, unsigned Total_Cells, int *Channel_Num_Cells, double *dx){
    
    kl=3;
    ku=3;
    ldab=2*kl+ku+1;
    
    Num_Channel_=Num_Channel;    
    Num_Reservoir_=Num_Reservoir;
    Total_Cells_=Total_Cells;
    
    ab=new double [ldab*2*Total_Cells_];
    a=new double [2*Num_Reservoir_*2*Num_Reservoir_];
    rhs_1=new double[2*Total_Cells_*3];
    rhs_2=new double [2*Num_Reservoir_];
    Reservoir_P=new double [Num_Reservoir_];
    Reservoir_Phi=new double [Num_Reservoir_];
    
    Channel_Num_Cells_=new int [Num_Channel_];
    dx_=new double [Num_Channel_];
    for(int i=0; i< Num_Channel_; i++)
    {
        Channel_Num_Cells_[i]=Channel_Num_Cells[i];
        dx_[i]=dx[i];
    }
    
    LHS_1=new double* [ldab];
    for(int i=0; i<ldab; i++)
   {
       LHS_1[i]=new double [2*Total_Cells_];
   }

   LHS_2=new double* [2*Num_Reservoir_];
   for(int i=0; i<2*Num_Reservoir_; i++)
   {
       LHS_2[i]=new double [2*Num_Reservoir_];
   }

   LHS_3=new double* [Num_Reservoir_];
   for(int i=0; i<Num_Reservoir_; i++)
   {
       LHS_3[i]=new double [Num_Reservoir_];
   }

   P=new double* [Total_Cells_];
   Phi=new double* [Total_Cells_];
   for(int i=0; i<Total_Cells_; i++)
   {
       P[i]=new double [3];
       Phi[i]=new double [3];
   }
  
   rhs_p_phi=new double[2*Total_Cells_];

}

Lapack::~Lapack(){
    
    delete [] Channel_Num_Cells_;
    delete [] dx_;
    delete [] ab;
    delete [] a;
    delete [] rhs_1;
    delete [] rhs_2;
    delete [] Reservoir_P;
    delete [] Reservoir_Phi;
       
    for(int i=0; i<ldab; i++)
    {
        delete [] LHS_1[i];
    }
    delete [] LHS_1;
        
    for(int i=0; i<2*Num_Reservoir_; i++)
    {
        delete [] LHS_2[i];
    }
    delete [] LHS_2;

     for(int i=0; i<Num_Reservoir_; i++)
    {
        delete [] LHS_3[i];
    }
    delete [] LHS_3;
    
    for(int i=0; i<Total_Cells_; i++)
    {
        delete [] P[i];
        delete [] Phi[i];
    }
    delete [] P;
    delete [] Phi;

    delete [] rhs_p_phi; 
    
}


void Lapack::dgesv(double **LHS_2, double *rhs_2, int n){
    
    size=n;
    nrhs = 1;
    lda = n;
    ldb = n;
    
    for (int j=0; j<lda; j++) //for column
    {
        for (int i=0; i<lda; i++) //for row
        {
            a[i+j*lda] = LHS_2[i][j];
        }
    }
        
    ipiv = new int[n]; // Allocate memory to store pivot information
          
    // Now the function call
    dgesv_(&size, &nrhs, a, &lda, ipiv, rhs_2, &ldb, &info);
    cout<<"message from dgesv_: info="<<info<<endl;
   
    delete [] ipiv;
    return;
}

void Lapack::dgbsv(double **LHS_1, double *rhs_1, int n, int nrhs_){
    
    size=n;
    nrhs=nrhs_;
    ldb=n;
    ipiv = new int[n];
    
    for (int j=0; j<n; j++) //for column
    {
        for (int i=0; i<ldab; i++) //for row
        {
            ab[i+j*ldab] = LHS_1[i][j];
        }
    }
    
    // Now the function call
    
    dgbsv_(&size, &kl, &ku, &nrhs, ab, &ldab, ipiv, rhs_1, &ldb, &info);
    cout<<"message from dgbsv_: info="<<info<<endl;
    
    delete [] ipiv;    
    return;
}

void Lapack::dgtsv(double *DL, double *D, double *DU, double *B, int n_, int nrhs_){
    
    size=n_;
    nrhs=nrhs_;
    ldb=n_;
    
    // Now the function call
    dgtsv_(&size, &nrhs, DL, D, DU, B, &ldb, &info);
    //cout<<"message from dgtsv_: info="<<info<<endl;

    return;
}

void Lapack::find_three_P_Phi(double *A1_x_, double *B1_x_, double *A2_x_, double *B2_x_, double *f_1_x_, double *f_2_x_, double *S_x_, double *RHS_,  double *dx, double *Q1, double *Q2, double *Q3, double *I1, double *I2, double *I3 ){
    
    ////setting to zero
    for(int j=0; j<2*Total_Cells_; j++)
    {
        for(int i=0; i<ldab; i++)
        {
            LHS_1[i][j]=0.;
        }

    rhs_1[j]=0.;
    rhs_1[2*Total_Cells_+j]=0.;
    rhs_1[4*Total_Cells_+j]=0.;
    
    }
    
    ////Preparing rhs_1:
    for(int i=0; i<2*Total_Cells_; i++)
    {
        rhs_1[i]=RHS_[i];
        rhs_1[2*Total_Cells_+i]=0.;
        rhs_1[4*Total_Cells_+i]=0.;
    }
    
    cell_counter=0;
    for(int i=0; i<Num_Channel_; i++)
    {
        //Pressure 
        rhs_1[cell_counter]=0.;
        rhs_1[2*Total_Cells_+cell_counter]=0.;
        rhs_1[4*Total_Cells_+cell_counter]=0.;
        
        //Potential
        rhs_1[cell_counter+1]=0.;
        rhs_1[2*Total_Cells_+cell_counter+1]=0.;
        rhs_1[4*Total_Cells_+cell_counter+1]=0.;
        
        cell_counter+=(Channel_Num_Cells_[i]+2)*2;
        
        //Pressure
        rhs_1[cell_counter-2]=0.;
        rhs_1[2*Total_Cells_+cell_counter-2]=1.;
        rhs_1[4*Total_Cells_+cell_counter-2]=0.;
        
        //Potential
        rhs_1[cell_counter-1]=0.;
        rhs_1[2*Total_Cells_+cell_counter-1]=0.;
        rhs_1[4*Total_Cells_+cell_counter-1]=1.; 
    }

    ///////////////////////
    
    ////Preparing LHS_1
    cell_counter=0;
    face_counter=0;
    
    for(int i=0; i<Num_Channel_; i++){
        
        for(int j=1; j<Channel_Num_Cells_[i]+1; j++){
            
            //super-diagonal 3
            LHS_1[kl][kl+2*cell_counter+2*j]=S_x_[face_counter+j]*B1_x_[face_counter+j]/(dx_[i]*dx_[i]);
            
            //super-diagonal 2
            LHS_1[kl+1][(kl-1)+2*cell_counter+2*j]=S_x_[face_counter+j]*A1_x_[face_counter+j]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+1][(kl-1)+2*cell_counter+(2*j+1)]=S_x_[face_counter+j]*B2_x_[face_counter+j]/(dx_[i]*dx_[i]);
           
            //super-diagonal 1
            LHS_1[kl+2][(kl-2)+2*cell_counter+2*j]=-S_x_[face_counter+j-1]*B1_x_[face_counter+j-1]/(dx_[i]*dx_[i]) -S_x_[face_counter+j]*B1_x_[face_counter+j]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+2][(kl-2)+2*cell_counter+(2*j+1)]=S_x_[face_counter+j]*A2_x_[face_counter+j]/(dx_[i]*dx_[i]);
            
            //main diagonal non-boundary elements
            LHS_1[kl+3][2*cell_counter+2*j]=-S_x_[face_counter+j-1]*A1_x_[face_counter+j-1]/(dx_[i]*dx_[i]) - S_x_[face_counter+j]*A1_x_[face_counter+j]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+3][2*cell_counter+(2*j+1)]=-S_x_[face_counter+j-1]*B2_x_[face_counter+j-1]/(dx_[i]*dx_[i]) -S_x_[face_counter+j]*B2_x_[face_counter+j]/(dx_[i]*dx_[i]);
            
            //sub-diagonal 1
            LHS_1[kl+4][2*cell_counter+2*j-1]=S_x_[face_counter+j-1]*B1_x_[face_counter+j-1]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+4][2*cell_counter+(2*j+1)-1]=-S_x_[face_counter+j-1]*A2_x_[face_counter+j-1]/(dx_[i]*dx_[i]) - S_x_[face_counter+j]*A2_x_[face_counter+j]/(dx_[i]*dx_[i]);
            
        
            //subdiagonal 2
            LHS_1[kl+5][2*cell_counter+2*j-2]= S_x_[face_counter+j-1]*A1_x_[face_counter+j-1]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+5][2*cell_counter+(2*j+1)-2]= S_x_[face_counter+j-1]*B2_x_[face_counter+j-1]/(dx_[i]*dx_[i]);
            
            //subdiagonal 3
            LHS_1[kl+6][2*cell_counter+2*j-2]= S_x_[face_counter+j-1]*A2_x_[face_counter+j-1]/(dx_[i]*dx_[i]);
        }
        
        //setting Boundary elements for main diagonal row
        LHS_1[kl+3][2*cell_counter]=1., LHS_1[kl+3][2*cell_counter+1]=1.;
        
        cell_counter+=Channel_Num_Cells_[i]+2;
        face_counter+=Channel_Num_Cells_[i]+1;
        
        LHS_1[kl+3][2*cell_counter-2]=1., LHS_1[kl+3][2*cell_counter-1]=1.;
    }
    
    //some particular elements to be set separately
    LHS_1[kl+4][0]=0.;
    
    ////Solving System
    dgbsv(LHS_1, rhs_1, 2*Total_Cells_, 3);
    
    ////Transferring solutions to P and Phi arrays
    for(int i=0; i<Total_Cells_; i++)
    {
        P[i][0]=rhs_1[2*i];
        P[i][1]=rhs_1[2*Total_Cells_+2*i];
        P[i][2]=rhs_1[4*Total_Cells_+2*i];
        
        Phi[i][0]=rhs_1[2*i+1];
        Phi[i][1]=rhs_1[2*Total_Cells_+2*i+1];
        Phi[i][2]=rhs_1[4*Total_Cells_+2*i+1];

    }
   
    ////Calculating Q and I
    cell_counter=0;
    face_counter=0;
    
    for(int i=0; i<Num_Channel_; i++)
    {
        ////Flow rate
        Q1[i]=A1_x_[face_counter+Channel_Num_Cells_[i]/2]* (P[cell_counter+Channel_Num_Cells_[i]/2+1][0]-P[cell_counter+Channel_Num_Cells_[i]/2][0])/dx_[i]+ B1_x_[face_counter+Channel_Num_Cells_[i]/2]* (Phi[cell_counter+Channel_Num_Cells_[i]/2+1][0]-Phi[cell_counter+Channel_Num_Cells_[i]/2][0])/dx_[i]-f_1_x_[face_counter+Channel_Num_Cells_[i]/2];
    	Q1[i] *= S_x_[face_counter+Channel_Num_Cells_[i]/2];
    
        Q2[i]=A1_x_[face_counter+Channel_Num_Cells_[i]/2]* (P[cell_counter+Channel_Num_Cells_[i]/2+1][1]-P[cell_counter+Channel_Num_Cells_[i]/2][1])/dx_[i]+ B1_x_[face_counter+Channel_Num_Cells_[i]/2]* (Phi[cell_counter+Channel_Num_Cells_[i]/2+1][1]-Phi[cell_counter+Channel_Num_Cells_[i]/2][1])/dx_[i];
	Q2[i] *= S_x_[face_counter+Channel_Num_Cells_[i]/2];
        
        Q3[i]=A1_x_[face_counter+Channel_Num_Cells_[i]/2]* (P[cell_counter+Channel_Num_Cells_[i]/2+1][2]-P[cell_counter+Channel_Num_Cells_[i]/2][2])/dx_[i]+ B1_x_[face_counter+Channel_Num_Cells_[i]/2]* (Phi[cell_counter+Channel_Num_Cells_[i]/2+1][2]-Phi[cell_counter+Channel_Num_Cells_[i]/2][2])/dx_[i];
	Q3[i] *= S_x_[face_counter+Channel_Num_Cells_[i]/2];
        
        ////Current
        I1[i]=A2_x_[face_counter+Channel_Num_Cells_[i]/2]*(P[cell_counter+Channel_Num_Cells_[i]/2+1][0]-P[cell_counter+Channel_Num_Cells_[i]/2][0])/dx_[i]+ B2_x_[face_counter+Channel_Num_Cells_[i]/2]* (Phi[cell_counter+Channel_Num_Cells_[i]/2+1][0]-Phi[cell_counter+Channel_Num_Cells_[i]/2][0])/dx_[i] -f_2_x_[face_counter+Channel_Num_Cells_[i]/2];
	I1[i] *= S_x_[face_counter+Channel_Num_Cells_[i]/2];
        
        I2[i]=A2_x_[face_counter+Channel_Num_Cells_[i]/2]* (P[cell_counter+Channel_Num_Cells_[i]/2+1][1]-P[cell_counter+Channel_Num_Cells_[i]/2][1])/dx_[i]+ B2_x_[face_counter+Channel_Num_Cells_[i]/2]* (Phi[cell_counter+Channel_Num_Cells_[i]/2+1][1]-Phi[cell_counter+Channel_Num_Cells_[i]/2][1])/dx_[i];
	I2[i] *= S_x_[face_counter+Channel_Num_Cells_[i]/2];
        
        I3[i]=A2_x_[face_counter+Channel_Num_Cells_[i]/2]* (P[cell_counter+Channel_Num_Cells_[i]/2+1][2]-P[cell_counter+Channel_Num_Cells_[i]/2][2])/dx_[i]+ B2_x_[face_counter+Channel_Num_Cells_[i]/2]* (Phi[cell_counter+Channel_Num_Cells_[i]/2+1][2]-Phi[cell_counter+Channel_Num_Cells_[i]/2][2])/dx_[i];
	I3[i] *= S_x_[face_counter+Channel_Num_Cells_[i]/2];
        
        cell_counter+=Channel_Num_Cells_[i]+2;
        face_counter+=Channel_Num_Cells_[i]+1;

        //cout<<i<<setw(20)<<Q1[i]<<setw(20)<<Q2[i]<<setw(20)<<Q3[i]<<endl;
        //cout<<i<<setw(20)<<I1[i]<<setw(20)<<I2[i]<<setw(20)<<I3[i]<<endl;
    }
    
    return;
}

void Lapack::find_reservoir_P_Phi(double *Q1, double *Q2, double *Q3, double *I1, double *I2, double *I3, int **Connectivity, int **Connecting_Channel, bool *Reservoir_Pressure_Type, bool *Reservoir_Potential_Type, double *Reservoir_Pressure_Value, double *Reservoir_Potential_Value){
    
    for(int i=0;i< 2*Num_Reservoir_; i++){
        
        for(int j=0;j< 2*Num_Reservoir_; j++){
            
            LHS_2[i][j]=0.;
        }
		
        rhs_2[i]=0.;
    }
    
    for(int i=0; i<Num_Reservoir_; i++){
        
        if(!Reservoir_Pressure_Type[i]){
            
            LHS_2[2*i][2*i]=1.;
            
            rhs_2[2*i]=Reservoir_Pressure_Value[i];
        }
        
        else
        {
            
            for(int j=0; j<Num_Reservoir_; j++)
            {
                //non-diagonal elements
                if(j!=i)
                {
                    LHS_2[2*i][2*j]= abs(Connectivity[i][j])*Q2[Connecting_Channel[i][j]];
                    
                    LHS_2[2*i][2*j+1]= abs(Connectivity[i][j])*Q3[Connecting_Channel[i][j]];
                }
                
                //diagonal elements
                
                LHS_2[2*i][2*i]+= - abs(Connectivity[i][j])*Q2[Connecting_Channel[i][j]];
                
                LHS_2[2*i][2*i+1]+= - abs(Connectivity[i][j])*Q3[Connecting_Channel[i][j]];
                
                rhs_2[2*i]+= - Connectivity[i][j]*Q1[Connecting_Channel[i][j]];                
            }
        }
        
        if(!Reservoir_Potential_Type[i])
        {
            LHS_2[2*i+1][2*i+1]=1.;
            
            rhs_2[2*i+1]=Reservoir_Potential_Value[i];
        }
        
        else
        {
            for(int j=0; j<Num_Reservoir_; j++)
            {
                
                //non-diagonal elements
                if(j!=i)
                {
                    LHS_2[2*i+1][2*j]= abs(Connectivity[i][j])*I2[Connecting_Channel[i][j]];
                    
                    LHS_2[2*i+1][2*j+1]= abs(Connectivity[i][j])*I3[Connecting_Channel[i][j]];
                }
                
                //diagonal elements
                LHS_2[2*i+1][2*i]+= - abs(Connectivity[i][j])*I2[Connecting_Channel[i][j]];
                
                LHS_2[2*i+1][2*i+1]+= - abs(Connectivity[i][j])*I3[Connecting_Channel[i][j]];
                
                rhs_2[2*i+1]+= - Connectivity[i][j]*I1[Connecting_Channel[i][j]];
            }
        }
    }
    
    dgesv(LHS_2, rhs_2, 2*Num_Reservoir_);
    
    for(int i=0; i<Num_Reservoir_; i++){
        
        Reservoir_P[i]=rhs_2[2*i];
        
        Reservoir_Phi[i]=rhs_2[2*i+1];
    }
    return;
}

void Lapack::find_final_P_Phi_Q_I(double *Q1, double *Q2, double *Q3, double *I1, double *I2, double *I3, int **Channel_End_Reservoir, double *Q, double *I, double *Pressure, double *Potential){
    
    for(int i=0; i<Num_Channel_; i++)
    {
        Q[i]= Q1[i] + (Reservoir_P[Channel_End_Reservoir[i][1]] - Reservoir_P[Channel_End_Reservoir[i][0]])*Q2[i] + (Reservoir_Phi[Channel_End_Reservoir[i][1]] - Reservoir_Phi[Channel_End_Reservoir[i][0]])*Q3[i];
        
        I[i]= I1[i] + (Reservoir_P[Channel_End_Reservoir[i][1]] - Reservoir_P[Channel_End_Reservoir[i][0]])*I2[i] + (Reservoir_Phi[Channel_End_Reservoir[i][1]] - Reservoir_Phi[Channel_End_Reservoir[i][0]])*I3[i];

        //cout<<i<<setw(20)<<Q[i]<<setw(20)<<I[i]<<endl;
    }
   
    cell_counter=0;
    
    for(int i=0; i<Num_Channel_; i++)
    {
        for(int j=0; j<Channel_Num_Cells_[i]+2; j++)
        {
            Pressure[cell_counter+j]=P[cell_counter+j][0] + Reservoir_P[Channel_End_Reservoir[i][0]] + (Reservoir_P[Channel_End_Reservoir[i][1]] - Reservoir_P[Channel_End_Reservoir[i][0]])* P[cell_counter+j][1] + (Reservoir_Phi[Channel_End_Reservoir[i][1]] - Reservoir_Phi[Channel_End_Reservoir[i][0]])*P[cell_counter+j][2];
        
            Potential[cell_counter+j]=Phi[cell_counter+j][0] + Reservoir_Phi[Channel_End_Reservoir[i][0]]+ (Reservoir_P[Channel_End_Reservoir[i][1]] - Reservoir_P[Channel_End_Reservoir[i][0]])* Phi[cell_counter+j][1] + (Reservoir_Phi[Channel_End_Reservoir[i][1]] - Reservoir_Phi[Channel_End_Reservoir[i][0]])*Phi[cell_counter+j][2];

            //cout.precision(16);
	    //cout<<j<<setw(20)<<Potential[cell_counter+j]<<endl;
        }
        
        cell_counter+=Channel_Num_Cells_[i]+2;
    }
    
    return;
}


void Lapack::P_Phi(double *A1_x_, double *B1_x_, double *A2_x_, double *B2_x_, double *RHS_, double *dx, int **Channel_End_Reservoir, double *f_1_x_, double *f_2_x_,  double *Q, double *I, double *Pressure, double *Potential){
 for(int j=0; j<2*Total_Cells_; j++)
    {
        for(int i=0; i<ldab; i++)
        {
            LHS_1[i][j]=0.;
        }
    }
 ////Preparing LHS_1
    cell_counter=0;
    
    for(int i=0; i<Num_Channel_; i++){
        
        for(int j=1; j<Channel_Num_Cells_[i]+1; j++){
            
            //super-diagonal 3
            LHS_1[kl][kl+2*cell_counter+2*j]=B1_x_[cell_counter+j]/(dx_[i]*dx_[i]);
            
            //super-diagonal 2
            LHS_1[kl+1][(kl-1)+2*cell_counter+2*j]=A1_x_[cell_counter+j]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+1][(kl-1)+2*cell_counter+(2*j+1)]=B2_x_[cell_counter+j]/(dx_[i]*dx_[i]);
            
            //super-diagonal 1
            LHS_1[kl+2][(kl-2)+2*cell_counter+2*j]=-B1_x_[cell_counter+j-1]/(dx_[i]*dx_[i]) -B1_x_[cell_counter+j]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+2][(kl-2)+2*cell_counter+(2*j+1)]=A2_x_[cell_counter+j]/(dx_[i]*dx_[i]);
            
            //main diagonal non-boundary elements
            LHS_1[kl+3][2*cell_counter+2*j]=-A1_x_[cell_counter+j-1]/(dx_[i]*dx_[i]) - A1_x_[cell_counter+j]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+3][2*cell_counter+(2*j+1)]=-B2_x_[cell_counter+j-1]/(dx_[i]*dx_[i]) -B2_x_[cell_counter+j]/(dx_[i]*dx_[i]);
            
            //sub-diagonal 1
            LHS_1[kl+4][2*cell_counter+2*j-1]=B1_x_[cell_counter+j-1]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+4][2*cell_counter+(2*j+1)-1]=-A2_x_[cell_counter+j-1]/(dx_[i]*dx_[i]) - A2_x_[cell_counter+j]/(dx_[i]*dx_[i]);
            
        
            //subdiagonal 2
            LHS_1[kl+5][2*cell_counter+2*j-2]= A1_x_[cell_counter+j-1]/(dx_[i]*dx_[i]);
            
            LHS_1[kl+5][2*cell_counter+(2*j+1)-2]= B2_x_[cell_counter+j-1]/(dx_[i]*dx_[i]);
            
            //subdiagonal 3
            LHS_1[kl+6][2*cell_counter+2*j-2]= A2_x_[cell_counter+j-1]/(dx_[i]*dx_[i]);
        }
        
        //setting Boundary elements for main diagonal row
        LHS_1[kl+3][2*cell_counter]=1., LHS_1[kl+3][2*cell_counter+1]=1.;
        
        cell_counter+=Channel_Num_Cells_[i]+2;
        
        LHS_1[kl+3][2*cell_counter-2]=1., LHS_1[kl+3][2*cell_counter-1]=1.;
    }
    
    //some particular elements to be set separately
    LHS_1[kl+4][0]=0.;
    
    ////Preparing rhs_1:
     for(int j=0; j<2*Total_Cells_; j++)
    {  
      rhs_p_phi[j]=0.;
    }

    for(int i=0; i<2*Total_Cells_; i++)
    {
        rhs_p_phi[i]=RHS_[i];
    }
  
    
    cell_counter=0;
    for(int i=0; i<Num_Channel_; i++)
    {
        //Pressure 
        rhs_p_phi[cell_counter]=Reservoir_P[Channel_End_Reservoir[i][0]];
       
        //Potential
        rhs_p_phi[cell_counter+1]=Reservoir_Phi[Channel_End_Reservoir[i][0]];
       
        
        cell_counter+=(Channel_Num_Cells_[i]+2)*2;
        
        //Pressure
        rhs_p_phi[cell_counter-2]=Reservoir_P[Channel_End_Reservoir[i][1]];
       
        
        //Potential
        rhs_p_phi[cell_counter-1]=Reservoir_Phi[Channel_End_Reservoir[i][1]];
    }

   ////Solving System
    dgbsv(LHS_1, rhs_p_phi, 2*Total_Cells_, 1);

    /*
    double p0, m;
    m=50/(1+dx_[0]);
    p0=50-m*dx_[0]/2;

    cout<<setprecision(16);
    
    for(int i=0; i<Total_Cells_; i++)
    {
       
      cout<<i<<setw(20)<<rhs_p_phi[2*i+1]<<setw(30)<<-m*(i-0.5)*dx_[0]+p0<<endl;
    }
    */
    //int tt;
    //cin>>tt;
    
    ////Transferring solutions to P and Phi arrays
    for(int i=0; i<Total_Cells_; i++)
    {
        Pressure[i]=rhs_p_phi[2*i];
        
        Potential[i]=rhs_p_phi[2*i+1];

        //cout<<setprecision(16);
        //cout<<i<<setw(20)<<Pressure[i]<<setw(30)<<Potential[i]<<endl;
    }
   

    for(int i=0;i<Num_Channel_;i++)
      {
	/////Flow rate
        Q[i]=A1_x_[cell_counter+Channel_Num_Cells_[i]/2]* (Pressure[cell_counter+Channel_Num_Cells_[i]/2+1]-Pressure[cell_counter+Channel_Num_Cells_[i]/2])/dx_[i]+ B1_x_[cell_counter+Channel_Num_Cells_[i]/2]* (Potential[cell_counter+Channel_Num_Cells_[i]/2+1]-Potential[cell_counter+Channel_Num_Cells_[i]/2])/dx_[i]-f_1_x_[cell_counter+Channel_Num_Cells_[i]/2];

	////Current
        I[i]=A2_x_[cell_counter+Channel_Num_Cells_[i]/2]*(Pressure[cell_counter+Channel_Num_Cells_[i]/2+1]-Pressure[cell_counter+Channel_Num_Cells_[i]/2])/dx_[i]+ B2_x_[cell_counter+Channel_Num_Cells_[i]/2]* (Potential[cell_counter+Channel_Num_Cells_[i]/2+1]-Potential[cell_counter+Channel_Num_Cells_[i]/2])/dx_[i] -f_2_x_[cell_counter+Channel_Num_Cells_[i]/2];
      }

    return;
}



void Lapack::solve(double *A1_x_, double *B1_x_, double *A2_x_, double *B2_x_, double *f_1_x_, double *f_2_x_, double *S_x_, double *RHS_, double *dx, int **Connectivity, int **Connecting_Channel, int **Channel_End_Reservoir, bool *Reservoir_Pressure_Type, bool *Reservoir_Potential_Type, double *Reservoir_Pressure_Value, double *Reservoir_Potential_Value, double *Q1, double *Q2, double *Q3, double *I1, double *I2, double *I3, double *Q, double *I, double *Pressure, double *Potential){
    
    find_three_P_Phi(A1_x_, B1_x_, A2_x_, B2_x_, f_1_x_, f_2_x_, S_x_, RHS_, dx, Q1, Q2, Q3, I1, I2, I3);
    
    find_reservoir_P_Phi(Q1, Q2, Q3, I1, I2, I3, Connectivity, Connecting_Channel,Reservoir_Pressure_Type, Reservoir_Potential_Type, Reservoir_Pressure_Value, Reservoir_Potential_Value);
    
    find_final_P_Phi_Q_I(Q1, Q2, Q3, I1, I2, I3, Channel_End_Reservoir, Q, I, Pressure, Potential);

    return;
}

void Lapack::find_reservoir_dc(double *Channel_Inlet_Flux_, double *Channel_Outlet_Flux_, double**three_inlet_flux_, double** three_outlet_flux_, bool*Reservoir_Pressure_Type, bool *Reservoir_Potential_Type, int **Connectivity, int **Connecting_Channel,double *Reservoir_Volume, double dt, double *reservoir_dc)
{

 for(int i=0;i< Num_Reservoir_; i++){
        
        for(int j=0;j< Num_Reservoir_; j++){
            
            LHS_3[i][j]=0.;
        }
		
        reservoir_dc[i]=0.;
    }
    
    for(int i=0; i<Num_Reservoir_; i++){
        
        if(!Reservoir_Pressure_Type[i] || !Reservoir_Potential_Type[i] ){
            
            LHS_3[i][i]=1.;
        }
        
        else
        {
	  LHS_3[i][i] = Reservoir_Volume[i]/dt;
            
	  for(int j=0; j<i; j++) //channel outlet flux which is entering the reservoir
            {

                //non-diagonal elements
                LHS_3[i][j]= -abs(Connectivity[j][i])*
		  three_outlet_flux_[Connecting_Channel[j][i]][1];
                
                //diagonal elements
                LHS_3[i][i]+= - abs(Connectivity[j][i])*
		  three_outlet_flux_[Connecting_Channel[j][i]][2];
                
		// rhs
                reservoir_dc[i]+= abs(Connectivity[j][i])*(Channel_Outlet_Flux_[Connecting_Channel[j][i]] + three_outlet_flux_[Connecting_Channel[j][i]][0]);                
            }

	  for(int j=i+1; j<Num_Reservoir_; j++)//channel inlet flux which is leaving the reservoir
	    {
	      //non-diagonal elements
	      LHS_3[i][j]= abs(Connectivity[i][j])*
		  three_inlet_flux_[Connecting_Channel[i][j]][2];
	      
	      //diagonal elements
	       LHS_3[i][i]+= abs(Connectivity[i][j])*
		  three_inlet_flux_[Connecting_Channel[i][j]][1];

	      //rhs
	       reservoir_dc[i] -= abs(Connectivity[i][j])*( Channel_Inlet_Flux_[Connecting_Channel[i][j]] + three_inlet_flux_[Connecting_Channel[i][j]][0]);     
	    }
  
        }//end of else
    }//end of for(i...)
          
    dgesv(LHS_3, reservoir_dc, Num_Reservoir_);
    
  return;
}


    
