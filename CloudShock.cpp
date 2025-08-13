#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "reconstruction.h"


void initilization();
void time_step();
void iteration();
void recon();
void eigen_value();
void flux();
void RK3();    
void LLF();
void bc();
void dimension_reconstruction();
void flux_reconstruction();
void hybrid();
void source();



double U[max_no_variable][points][points],F1[max_no_variable][points][points], flux_out[max_no_variable][points][points],flux_out2[max_no_variable][points][points],bx[points][points],u8[points][points],
       speed_sound[points][points],G1[max_no_variable][points][points];

double pressure[points][points],density[points][points],velocity_y[points][points],velocity_x[points][points],velocity[points][points],magnetic[points][points],velocity_z[points][points],
       ethalpy[points][points],x[points],U_new2[max_no_variable][points][points], Bz[points][points];

double delta_t=1e-04,delta_x,C_neg[points][points],C_pos[points][points],alpha,alpha_neg,afven_speed_x[points][points],afven_speed[points][points],temp1,slow_speed[points][points],fast_speed[points][points];


double lamda1[max_no_variable][points][points],lamda2[max_no_variable][points][points],U_1[max_no_variable][points][points],U_new[max_no_variable][points][points],
       U_pos[max_no_variable][points][points],U_neg[max_no_variable][points][points],C1_neg[points][points],C1_pos[points][points];

double lamda_max1[points][points],lamda_max2,min_lamda_3,delta_y,y[points],lamda_max4,afven_speed_y[points][points],slow_speed_y[points][points],fast_speed_y[points][points],lamda_max3[points][points];

double  NetdivBx[points][points], NetdivBy[points][points], divBy[points][points], divBx[points][points], Bp[points][points];

double UL1[max_no_variable][points][points],UL2[max_no_variable][points][points],UR1[max_no_variable][points][points],UR2[max_no_variable][points][points],S[max_no_variable][points][points];

double U_L[max_no_variable][points][points],U_R[max_no_variable][points][points],G_R[max_no_variable][points][points],G_L[max_no_variable][points][points], F_W[max_no_variable][points][points],G_W[max_no_variable][points][points];

double  W1=0.22,W2;

double lamda_max, CFL=.6, r,f,r_0=0.1,r_1=0.115;

int nj,ni,n1,n2,n3,n4,n5,n6,n7,n8;

double t=0, tf=0.06;

using namespace std;

int main()
{


FILE *density_file, *velocity_file, *pressure_file,*magnetic_file,*velocity_y_file,*magnetic3_file,*magneticz_file;
    printf("\n2-D SHOCK CLOUD INTERACTION WITH WENO-RUSANOV\n\n");


    density_file  = fopen("density.dat","w");
    velocity_file = fopen("velocity_x.dat","w");
    pressure_file = fopen("pressure.dat","w");
    magnetic_file = fopen("magnetic.dat","w");
  velocity_y_file = fopen("velocity_y.dat","w");
   magnetic3_file = fopen("magnetic_x.dat","w");
    magneticz_file = fopen("magnetic_p.dat","w"); 

	delta_x = L/cells;
    delta_y = L/cells;
    initilization();
	
   	   do{

        
	      RK3();
	     
	     // periodic();
	      time_step();
	    
		  //cout<<size<<endl;
	
          }
	   while(t<tf);
	   
	   
	   
	       for(i=8; i<=size-8; i=i+2)
	    {
	     for(j=8; j<=size-8; j=j+2)
	  {
		fprintf(density_file,"%lf\t%lf\t%lf\n",  x[j],y[i],density[i][j]);
		fprintf(velocity_file,"%lf\t%lf\t%lf\n", x[j],y[i],velocity_x[i][j]);
		fprintf(pressure_file,"%lf\t%lf\t%lf\n",  x[j],y[i],pressure[i][j]);
        fprintf(magnetic_file,"%lf\t%lf\t%lf\n",  x[j],y[i],magnetic[i][j]);
       	fprintf(velocity_y_file,"%lf\t%lf\t%lf\n", x[j],y[i],velocity_y[i][j]);
		fprintf(magnetic3_file,"%lf\t%lf\t%lf\n",  x[j],y[i],bx[i][j]);
        fprintf(magneticz_file,"%lf\t%lf\t%lf\n", x[j],y[i],Bp[i][j]);	   
	   	}
      	
    	fprintf(density_file,"\n");
		fprintf(velocity_file,"\n");
		fprintf(pressure_file,"\n");
        fprintf(magnetic_file,"\n");
        fprintf(magnetic3_file,"\n");
		fprintf(velocity_y_file,"\n");
	    fprintf(magneticz_file,"\n");
	
	
	
	
	}
    fclose(density_file);
    fclose(velocity_file);
    fclose(velocity_y_file);
	fclose(pressure_file);
    fclose(magnetic_file);  
    fclose(magnetic3_file);
     fclose(magneticz_file); 

cout<<"THE END OF THE PROGRAMME"<<endl;
}

void initilization()
{

 j=0;
   for(k=8; k<=size-8; k=k+2)                     //Here i for y axis, j for x-axis
    {
  
  
   y[k]= 0 + delta_y*(j);
  x[k]=  0+ delta_x*(j);
 
  j++;
  
}
  
  for(i=8; i<=size-8; i=i+2)                     //Here i for y axis, j for x-axis
    {
       
       
	   for(j=8; j<=size-8; j=j+2)
       {
	        r = sqrt((x[j]-0.8) * (x[j]-0.8) + (y[i]-0.5) * (y[i]-0.5));
	     
		    if(x[j]<0.6)
	        {
			
	     pressure[i][j] = 167.35;
	
	    density[i][j] = 3.87; 	         
    
	    magnetic[i][j] = 2.18; 
  	
	    velocity_x[i][j] = 0;
      
	    velocity_y[i][j] = 0;
      
	  	velocity_z[i][j] = 0;
	  
	  	Bz[i][j]         = -2.18 ;
	  
	    bx[i][j]        = 0 ;
	  
}
 

	  
	 else if(r<0.15)
	 {
	     pressure[i][j] = 1.0;
	
	    density[i][j] = 10.0; 	         
    
	    magnetic[i][j] = 0.56; 
  	
	    velocity_x[i][j] = -11.25;
      
	    velocity_y[i][j] = 0;
      
	  	velocity_z[i][j] = 0;
	  
	  	Bz[i][j]         = 0.56 ;
	  
	    bx[i][j]        = 0 ;
	 
	 
	 
	 }
		else 
	{
	
			
	     pressure[i][j] = 1;
	
	    density[i][j] = 1; 	         
    
	    magnetic[i][j] = 0.56; 
  	
	    velocity_x[i][j] = -11.25;
      
	    velocity_y[i][j] = 0;
      
	  	velocity_z[i][j] = 0;
	  
	  	Bz[i][j]         = 0.56 ;
	  
	    bx[i][j]        = 0 ;
	  
	 } 
	 }	         
}

	 ////////////////////////////////////////////////////////////////////////////////
	 for(i=8; i<=size-8; i=i+2)                                    //need to be checked 
    {

     for(j=8; j<=size-8; j=j+2)                                    //need to be checked 
    {
	
	

	
    speed_sound[i][j]    = sqrt(GAMMA*pressure[i][j]/density[i][j]);
	 
	ethalpy[i][j]    =     0.5*velocity_x[i][j]*velocity_x[i][j] +0.5*velocity_y[i][j]*velocity_y[i][j]+0.5*velocity_z[i][j]*velocity_z[i][j]+ speed_sound[i][j]*
                        speed_sound[i][j]/(GAMMA-1)+ 1.0*((Bz[i][j]*Bz[i][j] + magnetic[i][j]*magnetic[i][j]+bx[i][j]*bx[i][j])/density[i][j]);

    velocity[i][j]= sqrt(velocity_x[i][j]*velocity_x[i][j]+velocity_y[i][j]*velocity_y[i][j]+velocity_z[i][j]*velocity_z[i][j]);
    
	}
}
///////////////////////////////////////////////////////////////////////////////////ghost cells entry///////////////////////////////////////////////////////////////////////////
nj=size-6;
ni=size-6;

n1=size-10;
n2=size-10;

n3=size-4;
n4=size-2;
n5=size-12;
n6=size-14;
n7=size-16;
n8=size;

/////////////////////boundary cells extrapolation from interior....//////////
for(j=8; j<=size-8; j=j+2)
    {


  		pressure[ni][j]      =   pressure[size-8][j];
        density[ni][j]       =   density[size-8][j];
        ethalpy[ni][j]       =  ethalpy[size-8][j];
        velocity_x[ni][j]    =  velocity_x[size-8][j];
		velocity_y[ni][j]    =  velocity_y[size-8][j];
		velocity_z[ni][j]    = 	velocity_z[size-8][j];
          magnetic[ni][j]    =    magnetic[size-8][j];
          Bz[ni][j]          =   Bz[size-8][j];
          bx[ni][j]          =  bx[size-8][j] ;
	    speed_sound[ni][j]   =  speed_sound[size-8][j];           
    
	  }
     
    
 
    
   for(j=8; j<=size-8; j=j+2)
    {

 
	
			
		pressure[4][j]     =   pressure[8][j];
        density[4][j]     =   density[8][j];
        ethalpy[4][j]      =  ethalpy[8][j];
        velocity_x[4][j]    =  velocity_x[8][j]  ;
		velocity_y[4][j]    =  velocity_y[8][j] ;
		velocity_z[4][j]    = 	velocity_z[8][j] ;
          magnetic[4][j]    =    magnetic[8][j];
          Bz[4][j]          =    Bz[8][j]  ;
           bx[4][j]          =  bx[8][j]  ;
	    speed_sound[4][j]=     speed_sound[8][j];           
    //

    }
     
	
for(i=8; i<=size-8; i=i+2)
    {
    
        pressure[i][6]      =   pressure[i][8];
        density[i][6]       =   density[i][8];
        ethalpy[i][6]       =  ethalpy[i][8];
        velocity_x[i][6]    =  velocity_x[i][8];  
		velocity_y[i][6]    =  velocity_y[i][8] ;
		velocity_z[i][6]    = 	velocity_z[i][8]; 
          magnetic[i][6]    =    magnetic[i][8];
          Bz[i][6]          =    Bz[i][8];  
           bx[i][6]         =  bx[i][8] ; 
	    speed_sound[i][6]   =     speed_sound[i][8];           
          
	  
}

for(i=8; i<=size-8; i=i+2)
    {
  
        pressure[i][size-6]      =   pressure[i][size-8];
        density[i][size-6]       =   density[i][size-8];
        ethalpy[i][size-6]       =  ethalpy[i][size-8];
        velocity_x[i][size-6]    =  velocity_x[i][size-8];  
		velocity_y[i][size-6]    =  velocity_y[i][size-8] ;
		velocity_z[i][size-6]    = 	velocity_z[i][size-8]; 
          magnetic[i][size-6]    =   magnetic[i][size-8];
          Bz[i][size-6]          =  Bz[i][size-8];  
           bx[i][size-6]         =  bx[i][size-8] ; 
	    speed_sound[i][size-6]   =     speed_sound[i][size-8];           
  
	  }



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// for j-4  and size-4 cells



   for(j=8; j<=size-8; j=j+2)                    ///n3=size-2....//
    {

   		pressure[n3][j]     =   pressure[size-8][j];
        density[n3][j]     =   density[size-8][j];
        ethalpy[n3][j]      =  ethalpy[size-8][j];
        velocity_x[n3][j]    =  velocity_x[size-8][j];
		velocity_y[n3][j]    =  velocity_y[size-8][j];
		velocity_z[n3][j]    = 	velocity_z[size-8][j];
          magnetic[n3][j]    =    magnetic[size-8][j];
          Bz[n3][j]          =    Bz[size-8][j] ;
          bx[n3][j]          =  bx[size-8][j] ;
	    speed_sound[n3][j]=     speed_sound[size-8][j];           
    
	  }
     
    
 

   for(j=8; j<=size-8; j=j+4)
    {
	
		pressure[4][j]     =   pressure[8][j];
        density[4][j]     =   density[8][j];
        ethalpy[4][j]      =  ethalpy[8][j];
        velocity_x[4][j]    =  velocity_x[8][j]  ;
		velocity_y[4][j]    =  velocity_y[8][j] ;
		velocity_z[4][j]    = 	velocity_z[8][j] ;
          magnetic[4][j]    =    magnetic[8][j];
          Bz[4][j]          =    Bz[8][j]  ;
           bx[4][j]          =  bx[8][j]  ;
	    speed_sound[4][j]=     speed_sound[8][j];           
    //

    }
     
	
for(i=8; i<=size-8; i=i+2)
    {
    
   
		pressure[i][4]     =   pressure[i][8];
        density[i][4]     =   density[i][8];
        ethalpy[i][4]      =  ethalpy[i][8];
        velocity_x[i][4]    =  velocity_x[i][8];  
		velocity_y[i][4]    =  velocity_y[i][8] ;
		velocity_z[i][4]    = 	velocity_z[i][8]; 
          magnetic[i][4]    =    magnetic[i][8];
          Bz[i][4]          =    Bz[i][8];  
           bx[i][4]          =  bx[i][8] ; 
	    speed_sound[i][4]=     speed_sound[i][8];           
          
	  

}
for(i=8; i<=size-8; i=i+2)
    {
 


     
	 	pressure[i][n3]     =    pressure[i][nj];
        density[i][n3]      =    density[i][nj];
        ethalpy[i][n3]      =    ethalpy[i][nj];
        velocity_x[i][n3]   =   velocity_x[i][nj]  ;
		velocity_y[i][n3]   =   velocity_y[i][nj] ;
		velocity_z[i][n3]   = 	 velocity_z[i][nj] ;
          magnetic[i][n3]   =   magnetic[i][nj];
          Bz[i][n3]         =   Bz[i][nj]  ;
           bx[i][n3]        =  bx[i][nj]  ;
	   speed_sound[i][n3]   =  speed_sound[i][nj] ;          
    

}
            ////ok//
//////////////////////////////////////////////////////////////////////////////////////////////// size-2 and 2 cells.........////////////////////////////////
  for(j=8; j<=size-8; j=j+2)                    ///n3=size-2....//
    {

   		pressure[n4][j]     =   pressure[size-8][j];
        density[n4][j]     =   density[size-8][j];
        ethalpy[n4][j]      =  ethalpy[size-8][j];
        velocity_x[n4][j]    =  velocity_x[size-8][j];
		velocity_y[n4][j]    =  velocity_y[size-8][j];
		velocity_z[n4][j]    = 	velocity_z[size-8][j];
          magnetic[n4][j]    =    magnetic[size-8][j];
          Bz[n4][j]          =    Bz[size-8][j] ;
          bx[n4][j]          =  bx[size-8][j] ;
	    speed_sound[n4][j]=     speed_sound[size-8][j];           
    
	  }
     
    
 

   for(j=8; j<=size-8; j=j+4)
    {
	
		pressure[2][j]     =   pressure[8][j];
        density[2][j]     =   density[8][j];
        ethalpy[2][j]      =  ethalpy[8][j];
        velocity_x[2][j]    =  velocity_x[8][j]  ;
		velocity_y[2][j]    =  velocity_y[8][j] ;
		velocity_z[2][j]    = 	velocity_z[8][j] ;
          magnetic[2][j]    =    magnetic[8][j];
          Bz[2][j]          =    Bz[8][j]  ;
           bx[2][j]          =  bx[8][j]  ;
	    speed_sound[2][j]=     speed_sound[8][j];           
    //

    }
     
	
for(i=8; i<=size-8; i=i+2)
    {
    
   
		pressure[i][2]     =   pressure[i][8];
        density[i][2]     =   density[i][8];
        ethalpy[i][2]      =  ethalpy[i][8];
        velocity_x[i][2]    =  velocity_x[i][8];  
		velocity_y[i][2]    =  velocity_y[i][8] ;
		velocity_z[i][2]    = 	velocity_z[i][8]; 
          magnetic[i][2]    =    magnetic[i][8];
          Bz[i][2]          =    Bz[i][8];  
           bx[i][2]          =  bx[i][8] ; 
	    speed_sound[i][2]=     speed_sound[i][8];           
          
	  

}
for(i=8; i<=size-8; i=i+2)
    {
 


     
	 	pressure[i][n4]     =    pressure[i][nj];
        density[i][n4]      =    density[i][nj];
        ethalpy[i][n4]      =    ethalpy[i][nj];
        velocity_x[i][n4]   =   velocity_x[i][nj]  ;
		velocity_y[i][n4]   =   velocity_y[i][nj] ;
		velocity_z[i][n4]   = 	 velocity_z[i][nj] ;
          magnetic[i][n4]   =   magnetic[i][nj];
          Bz[i][n4]         =   Bz[i][nj]  ;
           bx[i][n4]        =  bx[i][nj]  ;
	   speed_sound[i][n4]   =  speed_sound[i][nj] ;          
    

}
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// for 0 and J cells /////////////////

    
   for(j=8; j<=size-8; j=j+2)
    {

          
		
		pressure[n8][j]     =   pressure[size-6][j];
        density[n8][j]     =   density[size-6][j];
        ethalpy[n8][j]      =  ethalpy[size-6][j];
        velocity_x[n8][j]    =  velocity_x[size-6][j];
		velocity_y[n8][j]    =  velocity_y[size-6][j];
		velocity_z[n8][j]    = 	velocity_z[size-6][j];
          magnetic[n8][j]    =    magnetic[size-6][j];
          Bz[n8][j]          =    Bz[size-6][j] ;
          bx[n8][j]          =  bx[size-6][j] ;
	    speed_sound[n8][j]=     speed_sound[size-6][j];           
    
	  }
    
    
for(j=8; j<=size-8; j=j+2)
    {

	
			
		pressure[0][j]     =   pressure[6][j];
        density[0][j]     =   density[6][j];
        ethalpy[0][j]      =  ethalpy[6][j];
        velocity_x[0][j]    =  velocity_x[6][j]  ;
		velocity_y[0][j]    =  velocity_y[6][j] ;
		velocity_z[0][j]    = 	velocity_z[6][j] ;
          magnetic[0][j]    =    magnetic[6][j];
          Bz[0][j]          =    Bz[6][j]  ;
           bx[0][j]          =  bx[6][j]  ;
	    speed_sound[0][j]=     speed_sound[6][j];           
    //

    }
    
	
for(i=8; i<=size-8; i=i+2)
    {
    
        
		
		
		pressure[i][0]     =   pressure[i][6];
        density[i][0]     =   density[i][6];
        ethalpy[i][0]      =  ethalpy[i][6];
        velocity_x[i][0]    =  velocity_x[i][6];  
		velocity_y[i][0]    =  velocity_y[i][6] ;
		velocity_z[i][0]    = 	velocity_z[i][6]; 
          magnetic[i][0]    =    magnetic[i][6];
          Bz[i][0]          =    Bz[i][6];  
           bx[i][0]          =  bx[i][6] ; 
	    speed_sound[i][0]=     speed_sound[i][6];           
          
	  }


for(i=8; i<=size-8; i=i+2)
    {
    
	 	pressure[i][n8]     =    pressure[i][nj];
        density[i][n8]      =    density[i][nj];
        ethalpy[i][n8]      =    ethalpy[i][nj];
        velocity_x[i][n8]   =   velocity_x[i][nj]  ;
		velocity_y[i][n8]   =   velocity_y[i][nj] ;
		velocity_z[i][n8]   = 	 velocity_z[i][nj] ;
          magnetic[i][n8]   =   magnetic[i][nj];
          Bz[i][n8]         =   Bz[i][nj]  ;
           bx[i][n8]        =  bx[i][nj]  ;
	   speed_sound[i][n8]   =  speed_sound[i][nj] ;          
    

}
/////////////////////////////////////////////////////// ghost cells entry done /////////////////////////////////////////
 






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  
for(i=6; i<=size-6; i=i+2)
    {
for(j=6; j<=size-6; j=j+2)
    {
          
        U[0][i][j]= density[i][j];
      
	    U[1][i][j]= density[i][j]*velocity_x[i][j];
      
	    U[2][i][j]= density[i][j]*velocity_y[i][j];
       
	    U[3][i][j]=density[i][j]*velocity_z[i][j];
      
	    U[4][i][j]= (pressure[i][j]/(GAMMA-1)) + (0.5*density[i][j]*(velocity_x[i][j]*velocity_x[i][j]+velocity_y[i][j]*velocity_y[i][j]+velocity_z[i][j]*velocity_z[i][j]))+(0.5*mu*(Bz[i][j]*Bz[i][j] + magnetic[i][j]*magnetic[i][j]+bx[i][j]*bx[i][j]));
      
	    U[5][i][j]= magnetic[i][j];
      
	    U[6][i][j]= Bz[i][j];

            U[7][i][j]= bx[i][j];

   // if(i==10)
//cout<<U[0][i][j]<<"\t"<<t<<"\n";
 }
}
return;
                                                                // have to include G flux and source terms 
}


void recon()
{


for(i=6; i<=size-6; i=i+2)
    {

for(j=6; j<=size-6; j=j+2)
   {
    afven_speed_x[i][j]      =  (bx[i][j])/(sqrt(density[i][j]));
	
 
 	afven_speed[i][j]     = (sqrt((bx[i][j]*bx[i][j])+(magnetic[i][j]*magnetic[i][j]+Bz[i][j]*Bz[i][j]))/sqrt(density[i][j]));
	

	

	temp1 = (speed_sound[i][j]*speed_sound[i][j])+(afven_speed[i][j]*afven_speed[i][j]);
	

	
	slow_speed[i][j] = sqrt(0.5*((temp1)-sqrt((temp1*temp1)-(4*speed_sound[i][j]*speed_sound[i][j]*afven_speed_x[i][j]*afven_speed_x[i][j]))));
	
	
	fast_speed[i][j] = sqrt(0.5*((temp1)+sqrt((temp1*temp1)-(4*speed_sound[i][j]*speed_sound[i][j]*afven_speed_x[i][j]*afven_speed_x[i][j]))));
	 
	 
   }
   
}
 //////////////////////////////////////////////// For Y terms
for(i=6; i<=size-6; i=i+2)
    {

for(j=6; j<=size-6; j=j+2)
   {
    afven_speed_y[i][j]      =  (magnetic[i][j])/(sqrt(density[i][j]));
	
 
	 
	 temp1 = (speed_sound[i][j]*speed_sound[i][j])+(afven_speed[i][j]*afven_speed[i][j]);
	

	
	slow_speed_y[i][j] = sqrt(0.5*((temp1)-sqrt((temp1*temp1)-(4*speed_sound[i][j]*speed_sound[i][j]*afven_speed_y[i][j]*afven_speed_y[i][j]))));
	
	
	fast_speed_y[i][j] = sqrt(0.5*((temp1)+sqrt((temp1*temp1)-(4*speed_sound[i][j]*speed_sound[i][j]*afven_speed_y[i][j]*afven_speed_y[i][j]))));
	 
   }
   
}
 
return;

 }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void eigen_value()
{
       
   lamda_max2=0;
   
 for(i=6; i<=size-6; i=i+2)
    {
for(j=6; j<=size-6; j=j+2)
   {	
    lamda1[3][i][j] = velocity_x[i][j];
	lamda1[2][i][j] = -slow_speed[i][j] + velocity_x[i][j];
	lamda1[4][i][j] = slow_speed[i][j] + velocity_x[i][j];
    lamda1[6][i][j] = velocity_x[i][j] + fast_speed[i][j];
	lamda1[0][i][j] =  velocity_x[i][j] - fast_speed[i][j];
    lamda1[5][i][j] = afven_speed_x[i][j] + velocity_x[i][j];
	lamda1[1][i][j] = -afven_speed_x[i][j]+ velocity_x[i][j];
    lamda1[7][i][j] = velocity_x[i][j];                                  //Powell term
	
//	if(lamda_max2<lamda1)
  
}
  
   }
     
  lamda_max4=0; 
  
  
   for(i=6; i<=size-6; i=i+2)
    {
	 for(j=6; j<=size-6; j=j+2)
   {
   
    lamda2[3][i][j] = velocity_y[i][j];
	lamda2[2][i][j] = -slow_speed_y[i][j] + velocity_y[i][j];
	lamda2[4][i][j] = slow_speed_y[i][j] + velocity_y[i][j];
    lamda2[6][i][j] = velocity_y[i][j] + fast_speed_y[i][j];
	lamda2[0][i][j] =  velocity_y[i][j] - fast_speed_y[i][j];
    lamda2[5][i][j] = afven_speed_y[i][j] + velocity_y[i][j];
	lamda2[1][i][j] = -afven_speed_y[i][j]+ velocity_y[i][j];
    lamda2[7][i][j] = velocity_y[i][j];                                  //Powell term
}


}

   
  
    for(i=6; i<=size-6; i=i+2)
    {

    for(j=6; j<=size-6; j=j+2)
    {

	
	  if(lamda_max2<fabs(lamda1[6][i][j]))  
          lamda_max2=fabs(lamda1[6][i][j]);
 
 
   	  
}

  }
  
 for(i=6; i<=size-6; i=i+2)
    {

    for(j=6; j<=size-6; j=j+2)
    {

	
	  if(lamda_max2<fabs(lamda1[0][i][j]))  
          lamda_max2=fabs(lamda1[0][i][j]);
 
 
   	  
}

  }
  for(i=6; i<=size-6; i=i+2)
    {

    for(j=6; j<=size-6; j=j+2)
    {

	
	  if(lamda_max4<fabs(lamda2[6][i][j]))   
          lamda_max4=fabs(lamda2[6][i][j]);
      
  

}

  }
  
 for(i=6; i<=size-6; i=i+2)
    {

    for(j=6; j<=size-6; j=j+2)
    {

	
	  if(lamda_max4<fabs(lamda2[0][i][j]))   
          lamda_max4=fabs(lamda2[0][i][j]);
      
  

}

  }
 

if(lamda_max2>lamda_max4)
  lamda_max=lamda_max2;
  
  else if(lamda_max2<lamda_max4)
  lamda_max=lamda_max4;
  
  else if(lamda_max2==lamda_max4)
     lamda_max=lamda_max4;
  	
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void LLF()
 {
           
   
   
	  for(i=6; i<=size-6; i=i+2)
       
	   {
      for(j=6; j<=size-6; j=j+2)
       
	   {
    
	   lamda_max1[i][j]=0;
	  
	   for(k=0; k<8; k++)
       {

      if(lamda_max1[i][j]<fabs(lamda1[k][i][j]))
         
		 lamda_max1[i][j]=fabs(lamda1[k][i][j]);
         
		 
      	 
		 
		  }

         }  
        }
        
		 for(i=6; i<=size-6; i=i+2)
       
	   {
      for(j=6; j<=size-6; j=j+2)
       
	   {
    
	   lamda_max3[i][j]=0;
	  
	   for(k=0; k<8; k++)
       {

         
		 
 if(lamda_max3[i][j]<fabs(lamda2[k][i][j]))
         
		 lamda_max3[i][j]=fabs(lamda2[k][i][j]);
      	 
		 
		  }

          }  
        }
        
        
  for(i=6; i<=size-8; i=i+2)
{

for(j=6; j<=size-8; j=j+2)
{
 
   
        if(lamda_max1[i][j]>lamda_max1[i][j+2])
          
C_pos[i][j+1] =lamda_max1[i][j];

    else
          
C_pos[i][j+1] =lamda_max1[i][j+2];

}
}
        
for(i=6; i<=size-8; i=i+2)
{

for(j=6; j<=size-8; j=j+2)
{

    

    if(lamda_max3[i+2][j]>lamda_max3[i][j])
          
C1_pos[i+1][j] =lamda_max3[i+2][j];

    else
          
C1_pos[i+1][j] =lamda_max3[i][j];

 //   if(i==50)
//cout<<C1_pos[i+1][j]<<"\t"<<t<<"\n";

  
  
   }
}


 return;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
 



void flux()
 
 {
 
 
 
  for(i=6; i<=size-8; i=i+2)
{
	                                            
  for(j=6; j<=size-8; j=j+2)
    {

      for(k=0;k<8;k++)
      {
     
       
       	

       UR1[k][i+1][j]= (U[k][i+2][j]- U[k][i][j])*C1_pos[i+1][j];
   
       UL1[k][i][j+1]= (U[k][i][j+2]- U[k][i][j])*C_pos[i][j+1];
   
  
  //if(i==50) 
  //cout<<UR1[1][i+1][j]<<"\t"<<t<<"\n";

       }
  
    }

}

for(i=6; i<=size-8; i=i+2)
	{

  for(j=6; j<=size-8; j=j+2)
    {
	
       F1[0][i][j+1] = (density[i][j]*velocity_x[i][j] + density[i][j+2]*velocity_x[i][j+2])*0.5 - .5*(UL1[0][i][j+1]);
       
	    
      
	    
		G1[0][i+1][j] = (density[i][j]*velocity_y[i][j] + density[i+2][j]*velocity_y[i+2][j])*0.5 - .5*(UR1[0][i+1][j]);
       
	    
	        
             
	 	F1[1][i][j+1] = (density[i][j+2]*velocity_x[i][j+2]*velocity_x[i][j+2] + density[i][j]*velocity_x[i][j]*velocity_x[i][j])*0.5 

                  +(.5*(pressure[i][j]+pressure[i][j+2])+.25*mu*(((magnetic[i][j]*magnetic[i][j])+(bx[i][j]*bx[i][j]))+(magnetic[i][j+2]*magnetic[i][j+2])+(bx[i][j+2]*bx[i][j+2])+(Bz[i][j]*Bz[i][j])+(Bz[i][j+2]*Bz[i][j+2]))) 
                  
                  -.5*(UL1[1][i][j+1])-((bx[i][j]*bx[i][j])+(bx[i][j+2]*bx[i][j+2]))/2;

	
	    G1[2][i+1][j] = (density[i][j]*velocity_y[i][j]*velocity_y[i][j] + density[i+2][j]*velocity_y[i+2][j]*velocity_y[i+2][j])*0.5 

                  +(.5*(pressure[i][j]+pressure[i+2][j])+.25*mu*(((magnetic[i][j]*magnetic[i][j])+(bx[i][j]*bx[i][j])+(Bz[i][j]*Bz[i][j]))+(magnetic[i+2][j]*magnetic[i+2][j])+(bx[i+2][j]*bx[i+2][j])+(Bz[i+2][j]*Bz[i+2][j]))) 
                  
                  -.5*(UR1[2][i+1][j])-((magnetic[i+2][j]*magnetic[i+2][j])+(magnetic[i][j]*magnetic[i][j]))/2;

	
	
 
	   	G1[1][i+1][j] = ((density[i+2][j]*velocity_x[i+2][j]*velocity_y[i+2][j])+(density[i][j]*velocity_x[i][j]*velocity_y[i][j]))*.5 
                               -((magnetic[i][j]*bx[i][j])+magnetic[i+2][j]*bx[i+2][j])*mu*0.5 - .5*(UR1[1][i+1][j]);


        F1[2][i][j+1] = ((density[i][j+2]*velocity_x[i][j+2]*velocity_y[i][j+2])+(density[i][j]*velocity_x[i][j]*velocity_y[i][j]))*.5 
                               -((magnetic[i][j]*bx[i][j])+magnetic[i][j+2]*bx[i][j+2])*mu*0.5 - .5*(UL1[2][i][j+1]);


  
         F1[3][i][j+1] = ((density[i][j+2]*velocity_x[i][j+2]*velocity_z[i][j+2])+(density[i][j]*velocity_x[i][j]*velocity_z[i][j]))*.5 
                               -((Bz[i][j]*bx[i][j])+(Bz[i][j+2]*bx[i][j+2]))*mu*.5-.5*(UL1[3][i][j+1]);

        
	    G1[3][i+1][j] = ((density[i+2][j]*velocity_y[i+2][j]*velocity_z[i+2][j])+(density[i][j]*velocity_y[i][j]*velocity_z[i][j]))*.5 
                               -((Bz[i][j]*magnetic[i][j])+(Bz[i+2][j]*magnetic[i+2][j]))*mu*.5-.5*(UR1[3][i+1][j]);


		 
        F1[5][i][j+1] =((magnetic[i][j]*velocity_x[i][j])+(magnetic[i][j+2]*velocity_x[i][j+2]))*.5 - ((bx[i][j]*velocity_y[i][j])+(bx[i][j+2]*velocity_y[i][j+2]))*.5-.5*(UL1[5][i][j+1]);  

        G1[7][i+1][j] =((bx[i][j]*velocity_y[i][j])+(bx[i+2][j]*velocity_y[i+2][j]))*.5 - ((magnetic[i][j]*velocity_x[i][j])+(magnetic[i+2][j]*velocity_x[i+2][j]))*.5-.5*(UR1[7][i+1][j]);  

	  
 

        F1[6][i][j+1] = ((Bz[i][j]*velocity_x[i][j])+(Bz[i][j+2]*velocity_x[i][j+2]))*.5 - ((bx[i][j]*velocity_z[i][j])+(bx[i][j+2]*velocity_z[i][j+2]))*.5-.5*(UL1[6][i][j+1]);

        G1[6][i+1][j] = ((Bz[i][i]*velocity_y[i][j])+(Bz[i+2][j]*velocity_y[i+2][j]))*.5 - ((magnetic[i][j]*velocity_z[i][j])+(magnetic[i+2][j]*velocity_z[i+2][j]))*.5-.5*(UR1[6][i+1][j]);
    


	  
	  
	    F1[4][i][j+1] = ((density[i][j]*velocity_x[i][j]*ethalpy[i][j]) + (density[i][j+2]*
                      velocity_x[i][j+2]*ethalpy[i][j+2]))*.5- (.5*mu*(bx[i][j]*(bx[i][j]*velocity_x[i][j]+(magnetic[i][j]*velocity_y[i][j])+(Bz[i][j]*velocity_z[i][j]))+(bx[i][j+2]*(bx[i][j+2]*velocity_x[i][j+2]+(magnetic[i][j+2]*velocity_y[i][j+2])+(Bz[i][j+2]*velocity_z[i][j+2])))))-.5*(UL1[4][i][j+1]); 

  
        G1[4][i+1][j] = ((density[i][j]*velocity_y[i][j]*ethalpy[i][j]) + (density[i+2][j]*velocity_y[i+2][j]*ethalpy[i+2][j]))*.5- (.5*mu*(magnetic[i][j]*(bx[i][j]*velocity_x[i][j]+(magnetic[i][j]*velocity_y[i][j])+(Bz[i][j]*velocity_z[i][j]))+(magnetic[i+2][j]*(bx[i+2][j]*velocity_x[i+2][j]+(magnetic[i+2][j]*velocity_y[i+2][j])+(Bz[i+2][j]*velocity_z[i+2][j])))))-.5*(UR1[4][i+1][j]); 

 
        
        F1[7][i][j+1] = -.5*(UL1[7][i][j+1]);


        G1[5][i+1][j] = -.5*(UR1[5][i+1][j]);
        
  
}

}
  
 

}

 
void RK3()
{

 //	bc();
 iteration();   
 //periodic();
				
    for(i=8; i<=size-12; i=i+2)
    {
 for(j=8; j<=size-12; j=j+2)
    {

    for(k=0; k<8; k++)
    {


   	U_1[k][i+2][j+2] = U[k][i+2][j+2] + (delta_t/delta_x)*(F1[k][i+2][j+1] - F1[k][i+2][j+3]) + (delta_t/delta_x)*(G1[k][i+1][j+2] - G1[k][i+3][j+2])- (delta_t)*(S[k][i+2][j+2]);

}



          density[i+2][j+2]     = U_1[0][i+2][j+2];
        velocity_x[i+2][j+2]    = U_1[1][i+2][j+2]/U_1[0][i+2][j+2];
        velocity_y[i+2][j+2]    = U_1[2][i+2][j+2]/U_1[0][i+2][j+2];
        velocity_z[i+2][j+2]    = U_1[3][i+2][j+2]/U_1[0][i+2][j+2];
          magnetic[i+2][j+2]    = U_1[5][i+2][j+2];
          Bz[i+2][j+2]          = U_1[6][i+2][j+2];   
           bx[i+2][j+2]          =   U_1[7][i+2][j+2];
	
	 pressure[i+2][j+2]    = (GAMMA-1)*(U_1[4][i+2][j+2] - 0.5*U_1[1][i+2][j+2]*U_1[1][i+2][j+2]/U_1[0][i+2][j+2]-0.5*U_1[2][i+2][j+2]*U_1[2][i+2][j+2]/U_1[0][i+2][j+2]-0.5*U_1[3][i+2][j+2]*U_1[3][i+2][j+2]/U_1[0][i+2][j+2]-0.5*(U_1[7][i+2][j+2]*U_1[7][i+2][j+2]+U_1[5][i+2][j+2]*U_1[5][i+2][j+2]+U_1[6][i+2][j+2]*U_1[6][i+2][j+2]));          
	
    speed_sound[i+2][j+2]    = sqrt(GAMMA*pressure[i+2][j+2]/density[i+2][j+2]);
	
	ethalpy[i+2][j+2]    =    0.5*velocity_x[i+2][j+2]*velocity_x[i+2][j+2] +0.5*velocity_y[i+2][j+2]*velocity_y[i+2][j+2]+0.5*velocity_z[i+2][j+2]*velocity_z[i+2][j+2]+ speed_sound[i+2][j+2]*
                         speed_sound[i+2][j+2]/(GAMMA-1)  + (1.0*(Bz[i+2][j+2]*Bz[i+2][j+2] + magnetic[i+2][j+2]*magnetic[i+2][j+2]+bx[i+2][j+2]*bx[i+2][j+2])/density[i+2][j+2]);






}
}



bc();
iteration();
 
   for(i=8; i<=size-12; i=i+2)
    {
 for(j=8; j<=size-12; j=j+2)
    {

    for(k=0; k<8; k++)
    {


   	U_new[k][i+2][j+2] = 0.75*U[k][i+2][j+2] +   0.25*U_1[k][i+2][j+2] + (0.25*delta_t/delta_x)*(F1[k][i+2][j+1] - F1[k][i+2][j+3]) + (0.25*delta_t/delta_x)*(G1[k][i+1][j+2] - G1[k][i+3][j+2]) - (0.25*delta_t)*(S[k][i+2][j+2]);

}



          density[i+2][j+2]     = U_new[0][i+2][j+2];
        velocity_x[i+2][j+2]    = U_new[1][i+2][j+2]/U_new[0][i+2][j+2];
        velocity_y[i+2][j+2]    = U_new[2][i+2][j+2]/U_new[0][i+2][j+2];
        velocity_z[i+2][j+2]    = U_new[3][i+2][j+2]/U_new[0][i+2][j+2];
          magnetic[i+2][j+2]    = U_new[5][i+2][j+2];
          Bz[i+2][j+2]          = U_new[6][i+2][j+2];   
           bx[i+2][j+2]          =   U_new[7][i+2][j+2];
	
	 pressure[i+2][j+2]    = (GAMMA-1)*(U_new[4][i+2][j+2] - 0.5*U_new[1][i+2][j+2]*U_new[1][i+2][j+2]/U_new[0][i+2][j+2]-0.5*U_new[2][i+2][j+2]*U_new[2][i+2][j+2]/U_new[0][i+2][j+2]-0.5*U_new[3][i+2][j+2]*U_new[3][i+2][j+2]/U_new[0][i+2][j+2]-0.5*(U_new[7][i+2][j+2]*U_new[7][i+2][j+2]+U_new[5][i+2][j+2]*U_new[5][i+2][j+2]+U_new[6][i+2][j+2]*U_new[6][i+2][j+2]));          
	
    speed_sound[i+2][j+2]    = sqrt(GAMMA*pressure[i+2][j+2]/density[i+2][j+2]);
	
	ethalpy[i+2][j+2]    =    0.5*velocity_x[i+2][j+2]*velocity_x[i+2][j+2] +0.5*velocity_y[i+2][j+2]*velocity_y[i+2][j+2]+0.5*velocity_z[i+2][j+2]*velocity_z[i+2][j+2]+ speed_sound[i+2][j+2]*
                         speed_sound[i+2][j+2]/(GAMMA-1)  + (1.0*(Bz[i+2][j+2]*Bz[i+2][j+2] + magnetic[i+2][j+2]*magnetic[i+2][j+2]+bx[i+2][j+2]*bx[i+2][j+2])/density[i+2][j+2]);






}
}  
 bc();
 iteration();
 
 
  for(i=8; i<=size-12; i=i+2)
    {
 for(j=8; j<=size-12; j=j+2)
    {

    for(k=0; k<8; k++)
    {


   	U_new2[k][i+2][j+2] = 0.33*U[k][i+2][j+2] +   0.67*U_new[k][i+2][j+2] + (0.67*delta_t/delta_x)*(F1[k][i+2][j+1] - F1[k][i+2][j+3]) + (0.67*delta_t/delta_x)*(G1[k][i+1][j+2] - G1[k][i+3][j+2]) -(0.67*delta_t)*(S[k][i+2][j+2]);
}



          density[i+2][j+2]     = U_new2[0][i+2][j+2];
        velocity_x[i+2][j+2]    = U_new2[1][i+2][j+2]/U_new2[0][i+2][j+2];
        velocity_y[i+2][j+2]    = U_new2[2][i+2][j+2]/U_new2[0][i+2][j+2];
        velocity_z[i+2][j+2]    = U_new2[3][i+2][j+2]/U_new2[0][i+2][j+2];
          magnetic[i+2][j+2]    = U_new2[5][i+2][j+2];
          Bz[i+2][j+2]          = U_new2[6][i+2][j+2];   
           bx[i+2][j+2]          =   U_new2[7][i+2][j+2];
	
	 pressure[i+2][j+2]    = (GAMMA-1)*(U_new2[4][i+2][j+2] - 0.5*U_new2[1][i+2][j+2]*U_new2[1][i+2][j+2]/U_new2[0][i+2][j+2]-0.5*U_new2[2][i+2][j+2]*U_new2[2][i+2][j+2]/U_new2[0][i+2][j+2]-0.5*U_new2[3][i+2][j+2]*U_new2[3][i+2][j+2]/U_new2[0][i+2][j+2]-0.5*(U_new2[7][i+2][j+2]*U_new2[7][i+2][j+2]+U_new2[5][i+2][j+2]*U_new2[5][i+2][j+2]+U_new2[6][i+2][j+2]*U_new2[6][i+2][j+2]));          
	
    speed_sound[i+2][j+2]    = sqrt(GAMMA*pressure[i+2][j+2]/density[i+2][j+2]);
	
	ethalpy[i+2][j+2]    =    0.5*velocity_x[i+2][j+2]*velocity_x[i+2][j+2] +0.5*velocity_y[i+2][j+2]*velocity_y[i+2][j+2]+0.5*velocity_z[i+2][j+2]*velocity_z[i+2][j+2]+ speed_sound[i+2][j+2]*
                         speed_sound[i+2][j+2]/(GAMMA-1)  + (1.0*(Bz[i+2][j+2]*Bz[i+2][j+2] + magnetic[i+2][j+2]*magnetic[i+2][j+2]+bx[i+2][j+2]*bx[i+2][j+2])/density[i+2][j+2]);

      Bp[i+2][j+2]= sqrt(bx[i+2][j+2]*bx[i+2][j+2]+magnetic[i+2][j+2]*magnetic[i+2][j+2]+ Bz[i+2][j+2]*Bz[i+2][j+2]) ;


if(i==50 && j==50)
cout<<"TIME\t"<<t<<"\tDIV.ERROR\t"<<abs(NetdivBx[i][j]+ NetdivBy[i][j])<<"\n";

 


}
}  
 
 
 for(i=10; i<=size-10; i=i+2)
    {
 for(j=10; j<=size-10; j=j+2)
    {

    for(k=0; k<8; k++)
    {
         U[k][i][j]=U_new2[k][i][j];  
 
 
}
}

}

 bc();		
}




void time_step()
{

 

	delta_t=CFL*delta_x/(lamda_max*2);


t=t+delta_t;


}
		
 void iteration()
	{
        
		recon();
    	eigen_value();
	    LLF();
        flux();
        dimension_reconstruction();
        flux_reconstruction();
        source();
		hybrid();  
      // periodic();   
		

	 }
		
		
void bc()
{


nj=size-6;
ni=size-6;

n1=size-10;
n2=size-10;

n3=size-4;
n4=size-2;
n5=size;
n6=size;

/////////////////////boundary cells extrapolation from interior....///////////////////////////////////////////////////


for(j=10; j<=size-10; j=j+2)
    {


  		pressure[8][j]     =   pressure[10][j];
        density[8][j]     =   density[10][j];
        ethalpy[8][j]      =  ethalpy[10][j];
        velocity_x[8][j]    =  velocity_x[10][j];
		velocity_y[8][j]    =  velocity_y[10][j];
		velocity_z[8][j]    = 	velocity_z[10][j];
          magnetic[8][j]    =   magnetic[10][j];
          Bz[8][j]          =    Bz[10][j] ;
          bx[8][j]          =  bx[10][j] ;
	    speed_sound[8][j]=     speed_sound[10][j];           
    
	  }
     
   


for(i=10; i<=size-10; i=i+2)
    {


  		pressure[i][8]     =   pressure[i][10];
        density[i][8]     =   density[i][10];
        ethalpy[i][8]      =  ethalpy[i][10];
        velocity_x[i][8]    =  velocity_x[i][10];
		velocity_y[i][8]    =  velocity_y[i][10];
		velocity_z[i][8]    = 	velocity_z[i][10];
          magnetic[i][8]    =    magnetic[i][10];
          Bz[i][8]          =   Bz[i][10] ;
          bx[i][8]          =  bx[i][10] ;
	    speed_sound[i][8]=     speed_sound[i][10];           
    
	  }
     


/////////////////////////////////////////////////////////////////// for Ghost cells next to the boundary............///////
 
	  for(j=10; j<=size-10; j=j+2)
    {                                                      //need to check on 31st august


  		pressure[size-8][j]     =   pressure[size-10][j];
        density[size-8][j]     =   density[size-10][j];
        ethalpy[size-8][j]      =  ethalpy[size-10][j];
        velocity_x[size-8][j]    =  velocity_x[size-10][j];
		velocity_y[size-8][j]    =  velocity_y[size-10][j];
		velocity_z[size-8][j]    = 	velocity_z[size-10][j];
          magnetic[size-8][j]    =    magnetic[size-10][j];
          Bz[size-8][j]          =    Bz[size-10][j] ;
          bx[size-8][j]          =  bx[size-10][j] ;
	    speed_sound[size-8][j]=     speed_sound[size-10][j];           
    
	  }
     


   // wrong right boundary..// maybe wright

for(i=10; i<=size-10; i=i+2)
    {

                      
  		pressure[i][size-10]     =   pressure[i][size-6];
        density[i][size-10]     =   density[i][size-6];
        ethalpy[i][size-10]      =  ethalpy[i][size-6];
        velocity_x[i][size-10]    =  velocity_x[i][size-6];
		velocity_y[i][size-10]    =  velocity_y[i][size-6];
		velocity_z[i][size-10]    = 	velocity_z[i][size-6];
          magnetic[i][size-10]    =   magnetic[i][size-6];
          Bz[i][size-10]          =  Bz[i][size-6] ;
          bx[i][size-10]          =  bx[i][size-6] ;
	    speed_sound[i][size-10]=     speed_sound[i][size-6];           
    
	  }
    
	  
///////////////////////////////////////////////////////////// ghost cells entry ///////////////// bc...........///////////
	  
   for(j=10; j<=size-10; j=j+2)
    {


  		pressure[ni][j]      =   pressure[size-8][j];
        density[ni][j]       =   density[size-8][j];
        ethalpy[ni][j]       =  ethalpy[size-8][j];
        velocity_x[ni][j]    =  velocity_x[size-8][j];
		velocity_y[ni][j]    =  velocity_y[size-8][j];
		velocity_z[ni][j]    = 	velocity_z[size-8][j];
          magnetic[ni][j]    =    magnetic[size-8][j];
          Bz[ni][j]          =    Bz[size-8][j] ;
          bx[ni][j]          =  bx[size-8][j] ;
	    speed_sound[ni][j]   =  speed_sound[size-8][j];           
    
	  }
     
    
 
    
   for(j=10; j<=size-10; j=j+2)
    {

 
	
			
		pressure[6][j]     =   pressure[8][j];
        density[6][j]     =   density[8][j];
        ethalpy[6][j]      =  ethalpy[8][j];
        velocity_x[6][j]    =  velocity_x[8][j]  ;
		velocity_y[6][j]    =  velocity_y[8][j] ;
		velocity_z[6][j]    = 	velocity_z[8][j] ;
          magnetic[6][j]    =    magnetic[8][j];
          Bz[6][j]          =    Bz[8][j]  ;
           bx[6][j]          =  bx[8][j]  ;
	    speed_sound[6][j]=     speed_sound[8][j];           
    //

    }
     
	
for(i=10; i<=size-10; i=i+2)
    {
    
        pressure[i][6]      =   pressure[i][8];
        density[i][6]       =   density[i][8];
        ethalpy[i][6]       =  ethalpy[i][8];
        velocity_x[i][6]    =  velocity_x[i][8];  
		velocity_y[i][6]    =  velocity_y[i][8] ;
		velocity_z[i][6]    = 	velocity_z[i][8]; 
          magnetic[i][6]    =    magnetic[i][8];
          Bz[i][6]          =    Bz[i][8];  
           bx[i][6]         =  bx[i][8] ; 
	    speed_sound[i][6]   =     speed_sound[i][8];           
          
	  
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// for j-2  and size-2 cells



   for(j=10; j<=size-10; j=j+2)                    ///n3=size-4....//cells.......//
    {

   		pressure[n3][j]     =   pressure[size-6][j];
        density[n3][j]     =   density[size-6][j];
        ethalpy[n3][j]      =  ethalpy[size-6][j];
        velocity_x[n3][j]    =  velocity_x[size-6][j];
		velocity_y[n3][j]    =  velocity_y[size-6][j];
		velocity_z[n3][j]    = 	velocity_z[size-6][j];
          magnetic[n3][j]    =    magnetic[size-6][j];
          Bz[n3][j]          =    Bz[size-6][j] ;
          bx[n3][j]          =  bx[size-6][j] ;
	    speed_sound[n3][j]=     speed_sound[size-6][j];           
    
	  }
     
    
 

   for(j=10; j<=size-10; j=j+2)
    {
	
		pressure[4][j]     =   pressure[6][j];
        density[4][j]     =   density[6][j];
        ethalpy[4][j]      =  ethalpy[6][j];
        velocity_x[4][j]    =  velocity_x[6][j]  ;
		velocity_y[4][j]    =  velocity_y[6][j] ;
		velocity_z[4][j]    = 	velocity_z[6][j] ;
          magnetic[4][j]    =    magnetic[6][j];
          Bz[4][j]          =    Bz[6][j]  ;
           bx[4][j]          =  bx[6][j]  ;
	    speed_sound[4][j]=     speed_sound[6][j];           
    //

    }
     
	
for(i=10; i<=size-10; i=i+2)
    {
    
   
		pressure[i][4]     =   pressure[i][6];
        density[i][4]     =   density[i][6];
        ethalpy[i][4]      =  ethalpy[i][6];
        velocity_x[i][4]    =  velocity_x[i][6];  
		velocity_y[i][4]    =  velocity_y[i][6] ;
		velocity_z[i][4]    = 	velocity_z[i][6]; 
          magnetic[i][4]    =    magnetic[i][6];
          Bz[i][4]          =    Bz[i][6];  
           bx[i][4]          =  bx[i][6] ; 
	    speed_sound[i][4]=     speed_sound[i][6];           
          
	  

}



            ////ok//

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// for 2 and size-2 cells /////////////////
  for(j=10; j<=size-10; j=j+2)                    ///n3=size-4....//cells.......//
    {

   		pressure[n4][j]     =   pressure[size-6][j];
        density[n4][j]     =   density[size-6][j];
        ethalpy[n4][j]      =  ethalpy[size-6][j];
        velocity_x[n4][j]    =  velocity_x[size-6][j];
		velocity_y[n4][j]    =  velocity_y[size-6][j];
		velocity_z[n4][j]    = 	velocity_z[size-6][j];
          magnetic[n4][j]    =    magnetic[size-6][j];
          Bz[n4][j]          =    Bz[size-6][j] ;
          bx[n4][j]          =  bx[size-6][j] ;
	    speed_sound[n4][j]=     speed_sound[size-6][j];           
    
	  }
     
    
 

   for(j=10; j<=size-10; j=j+2)
    {
	
		pressure[2][j]     =   pressure[6][j];
        density[2][j]     =   density[6][j];
        ethalpy[2][j]      =  ethalpy[6][j];
        velocity_x[2][j]    =  velocity_x[6][j]  ;
		velocity_y[2][j]    =  velocity_y[6][j] ;
		velocity_z[2][j]    = 	velocity_z[6][j] ;
          magnetic[2][j]    =    magnetic[6][j];
          Bz[2][j]          =    Bz[6][j]  ;
           bx[2][j]          =  bx[6][j]  ;
	    speed_sound[2][j]=     speed_sound[6][j];           
    //

    }
     
	
for(i=10; i<=size-10; i=i+2)
    {
    
   
		pressure[i][2]     =   pressure[i][6];
        density[i][2]     =   density[i][6];
        ethalpy[i][2]      =  ethalpy[i][6];
        velocity_x[i][2]    =  velocity_x[i][6];  
		velocity_y[i][2]    =  velocity_y[i][6] ;
		velocity_z[i][2]    = 	velocity_z[i][6]; 
          magnetic[i][2]    =    magnetic[i][6];
          Bz[i][2]          =    Bz[i][6];  
           bx[i][2]          =  bx[i][6] ; 
	    speed_sound[i][2]=     speed_sound[i][6];           
          
	  

}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
   for(j=10; j<=size-10; j=j+2)
    {

          
		
		pressure[n5][j]     =   pressure[size-6][j];
        density[n5][j]     =   density[size-6][j];
        ethalpy[n5][j]      =  ethalpy[size-6][j];
        velocity_x[n5][j]    =  velocity_x[size-6][j];
		velocity_y[n5][j]    =  velocity_y[size-6][j];
		velocity_z[n5][j]    = 	velocity_z[size-6][j];
          magnetic[n5][j]    =    magnetic[size-6][j];
          Bz[n5][j]          =    Bz[size-6][j] ;
          bx[n5][j]          =  bx[size-6][j] ;
	    speed_sound[n5][j]=     speed_sound[size-6][j];           
    
	  }
    
    
for(j=10; j<=size-10; j=j+2)
    {

	
			
		pressure[0][j]     =   pressure[6][j];
        density[0][j]     =   density[6][j];
        ethalpy[0][j]      =  ethalpy[6][j];
        velocity_x[0][j]    =  velocity_x[6][j]  ;
		velocity_y[0][j]    =  velocity_y[6][j] ;
		velocity_z[0][j]    = 	velocity_z[6][j] ;
          magnetic[0][j]    =    magnetic[6][j];
          Bz[0][j]          =    Bz[6][j]  ;
           bx[0][j]          =  bx[6][j]  ;
	    speed_sound[0][j]=     speed_sound[6][j];           
    //

    }
    
	
for(i=10; i<=size-10; i=i+2)
    {
    
        
		
		
		pressure[i][0]     =   pressure[i][6];
        density[i][0]     =   density[i][6];
        ethalpy[i][0]      =  ethalpy[i][6];
        velocity_x[i][0]    =  velocity_x[i][6];  
		velocity_y[i][0]    =  velocity_y[i][6] ;
		velocity_z[i][0]    = 	velocity_z[i][6]; 
          magnetic[i][0]    =    magnetic[i][6];
          Bz[i][0]          =    Bz[i][6];  
           bx[i][0]          =  bx[i][6] ; 
	    speed_sound[i][0]=     speed_sound[i][6];           
          
	  }






     for(i=8; i<=size-8; i=i+2)
{

for(k=0; k<8; k++)
    {

  
  		U[k][i][size-10] = U[k][i][size-8];
}
}

	 for(i=8; i<=size-8; i=i+2)
{

for(k=0; k<8; k++)
    {

  
  		U[k][i][8] = U[k][i][10];
}
}

     for(j=8; j<=size-8; j=j+2)
{

for(k=0; k<8; k++)
    {

  
  		U[k][size-8][j] = U[k][size-10][j];
}
}

 for(j=8; j<=size-8; j=j+2)
{

for(k=0; k<8; k++)
    {

  
  		U[k][8][j] = U[k][10][j];
}
}


return;

   // if(i==10)
//cout<<U[0][i][j]<<"\t"<<t<<"\n";
 }

	
void dimension_reconstruction()
{

WENO1(density);
WENO2(pressure);
WENO3(velocity_y);
WENO5(velocity_x);
WENO4(ethalpy);
WENO6(magnetic);
WENO7(bx);
WENO8(velocity_z);
WENO9(Bz);


}

void flux_reconstruction()
{


for(i=6; i<=size-8; i=i+2)
    {
for(j=6; j<=size-8; j=j+2)
    {
          
        U_L[0][i][j+1]= density_L[i][j+1];
      
	    U_L[1][i][j+1]= density_L[i][j+1]*velocity_xL[i][j+1];
      
	    U_L[2][i][j+1]= density_L[i][j+1]*velocity_yL[i][j+1];
       
	    U_L[3][i][j+1]=density_L[i][j+1]*velocity_zL[i][j+1];
      
	    U_L[4][i][j+1]= (pressure_L[i][j+1]/(GAMMA-1)) + (0.5*density_L[i][j+1]*(velocity_xL[i][j+1]*velocity_xL[i][j+1]+velocity_yL[i][j+1]*velocity_yL[i][j+1]+velocity_zL[i][j+1]*velocity_zL[i][j+1]))+(0.5*mu*(bz_L[i][j+1]*bz_L[i][j+1] + magnetic_L[i][j+1]*magnetic_L[i][j+1]+bx_L[i][j+1]*bx_L[i][j+1]));
      
	    U_L[5][i][j+1]= magnetic_L[i][j+1];
      
	    U_L[6][i][j+1]= bz_L[i][j+1];

            U_L[7][i][j+1]= bx_L[i][j+1];




           U_R[0][i][j-1]= density_R[i][j-1];
      
	    U_R[1][i][j-1]= density_R[i][j-1]*velocity_xR[i][j-1];
      
	    U_R[2][i][j-1]= density_R[i][j-1]*velocity_yR[i][j-1];
       
	    U_R[3][i][j-1]=density_R[i][j-1]*velocity_zR[i][j-1];
      
	    U_R[4][i][j-1]= (pressure_R[i][j-1]/(GAMMA-1)) + (0.5*density_R[i][j-1]*(velocity_xR[i][j-1]*velocity_xR[i][j-1]+velocity_yR[i][j-1]*velocity_yR[i][j-1]+velocity_zR[i][j-1]*velocity_zR[i][j-1]))+(0.5*mu*(bz_R[i][j-1]*bz_R[i][j-1] + magnetic_R[i][j-1]*magnetic_R[i][j-1]+bx_R[i][j-1]*bx_R[i][j-1]));
      
	    U_R[5][i][j-1]= magnetic_R[i][j-1];
      
	    U_R[6][i][j-1]= bz_R[i][j-1];

        U_R[7][i][j-1]= bx_R[i][j-1];


////////////////////////////////////////////////////////////////////// y-direction state variables //////////////////////////

        G_L[0][i+1][j]= d_L[i+1][j];
      
	    G_L[1][i+1][j]= d_L[i+1][j]*v_xL[i+1][j];
      
	    G_L[2][i+1][j]= d_L[i+1][j]*v_yL[i+1][j];
       
	    G_L[3][i+1][j]= d_L[i+1][j]*v_zL[i+1][j];
      
	    G_L[4][i+1][j]= (p_L[i+1][j]/(GAMMA-1)) + (0.5*d_L[i+1][j]*(v_xL[i+1][j]*v_xL[i+1][j]+v_yL[i+1][j]*v_yL[i+1][j]+v_zL[i+1][j]*v_zL[i+1][j]))+(0.5*mu*(BZ_L[i+1][j]*BZ_L[i+1][j] + BY_L[i+1][j]*BY_L[i+1][j]+BX_L[i+1][j]*BX_L[i+1][j]));
      
	    G_L[5][i+1][j]= BY_L[i+1][j];
      
	    G_L[6][i+1][j]= BZ_L[i+1][j];

        G_L[7][i+1][j]= BX_L[i+1][j];




        G_R[0][i-1][j]= d_R[i-1][j];
      
	    G_R[1][i-1][j]= d_R[i-1][j]*v_xR[i-1][j];
      
	    G_R[2][i-1][j]= d_R[i-1][j]*v_yR[i-1][j];
       
	    G_R[3][i-1][j]=  d_R[i-1][j]*v_zR[i-1][j];
      
	    G_R[4][i-1][j]= (p_R[i-1][j]/(GAMMA-1)) + (0.5*d_R[i-1][j]*(v_xR[i-1][j]*v_xR[i-1][j]+v_yR[i-1][j]*v_yR[i-1][j]+v_zR[i-1][j]*v_zR[i-1][j]))+(0.5*mu*(BZ_R[i-1][j]*BZ_R[i-1][j] + BY_R[i-1][j]*BY_R[i-1][j]+BX_R[i-1][j]*BX_R[i-1][j]));
      
	    G_R[5][i-1][j]= BY_R[i-1][j];
      
	    G_R[6][i-1][j]= BZ_R[i-1][j];

        G_R[7][i-1][j]= BX_R[i-1][j];



 }
}


for(i=6; i<=size-8; i=i+2)
    {
	                                            
  for(j=6; j<=size-8; j=j+2)
    {

      for(k=0;k<8;k++)
      {
     
       
       	

       UR2[k][i+1][j]= (G_R[k][i+1][j]- G_L[k][i+1][j])*C1_pos[i+1][j];
   
  
   
      

       UL2[k][i][j+1]= (U_R[k][i][j+1]- U_L[k][i][j+1])*C_pos[i][j+1];
   
  
  //if(i==50) 
  //cout<<UR1[1][i+1][j]<<"\t"<<t<<"\n";

}
  
}

}

for(i=6; i<=size-8; i=i+2)
	{

  for(j=6; j<=size-8; j=j+2)
    {
	
       F_W[0][i][j+1] = (density_R[i][j+1]*velocity_xR[i][j+1] + density_L[i][j+1]*velocity_xL[i][j+1])*0.5 - .5*(UL2[0][i][j+1]);
       
	    
      
	    
		G_W[0][i+1][j] = (d_L[i+1][j]*v_yL[i+1][j] + d_R[i+1][j]*v_yR[i+1][j])*0.5 - .5*(UR2[0][i+1][j]);
       
	    
	        
             
	 	F_W[1][i][j+1] = (density_R[i][j+1]*velocity_xR[i][j+1]*velocity_xR[i][j+1] + density_L[i][j+1]*velocity_xL[i][j+1]*velocity_xL[i][j+1])*0.5 

                  +(.5*(pressure_R[i][j+1]+pressure_L[i][j+1])+.25*mu*(((magnetic_L[i][j+1]*magnetic_L[i][j+1])+(bx_L[i][j+1]*bx_L[i][j+1]))+(magnetic_R[i][j+1]*magnetic_R[i][j+1])+(bx_R[i][j+1]*bx_R[i][j+1])+(bz_L[i][j+1]*bz_L[i][j+1])+(bz_R[i][j+1]*bz_R[i][j+1]))) 
                  
                  -.5*(UL2[1][i][j+1])-((bx_L[i][j+1]*bx_L[i][j+1])+(bx_R[i][j+1]*bx_R[i][j+1]))/2;

	
	    G_W[2][i+1][j] = (d_L[i+1][j]*v_yL[i+1][j]*v_yL[i+1][j] + d_R[i+1][j]*v_yR[i+1][j]*v_yR[i+1][j])*0.5 

                  +(.5*(p_L[i+1][j]+p_R[i+1][j])+.25*mu*(((BY_L[i+1][j]*BY_L[i+1][j])+(BX_L[i+1][j]*BX_L[i+1][j])+(BZ_R[i+1][j]*BZ_R[i+1][j]))+(BY_R[i+1][j]*BY_R[i+1][j])+(BX_R[i+1][j]*BX_R[i+1][j])+(BZ_L[i+1][j]*BZ_L[i+1][j]))) 
                  
                  -.5*(UR2[2][i+1][j])-((BY_L[i+1][j]*BY_L[i+1][j])+(BY_R[i+1][j]*BY_R[i+1][j]))/2;

	
	
 
	   	G_W[1][i+1][j] = ((d_L[i+1][j]*v_xL[i+1][j]*v_yL[i+1][j])+(d_R[i+1][j]*v_xR[i+1][j]*v_yR[i+1][j]))*.5 
                               -((BY_L[i+1][j]*BX_L[i+1][j])+BY_R[i+1][j]*BX_R[i+1][j])*mu*0.5 - .5*(UR2[1][i+1][j]);


        F_W[2][i][j+1] = ((density_L[i][j+1]*velocity_xL[i][j+1]*velocity_yL[i][j+1])+(density_R[i][j+1]*velocity_xR[i][j+1]*velocity_yR[i][j+1]))*.5 
                               -((magnetic_L[i][j+1]*bx_L[i][j+1])+magnetic_R[i][j+1]*bx_R[i][j+1])*mu*0.5 - .5*(UL2[2][i][j+1]);


  
         F_W[3][i][j+1] = ((density_L[i][j+1]*velocity_xL[i][j+1]*velocity_zL[i][j+1])+(density_R[i][j+1]*velocity_xR[i][j+1]*velocity_zR[i][j+1]))*.5 
                               -((bz_L[i][j+1]*bx_L[i][j+1])+(bz_R[i][j+1]*bx_R[i][j+1]))*mu*.5-.5*(UL2[3][i][j+1]);

        
	    G_W[3][i+1][j] = ((d_L[i+1][j]*v_yL[i+1][j]*v_zL[i+1][j])+(d_R[i+1][j]*v_yR[i+1][j]*v_zR[i+1][j]))*.5 
                               -((BZ_L[i+1][j]*BY_L[i+1][j])+(BZ_R[i+1][j]*BY_R[i+1][j]))*mu*.5-.5*(UR2[3][i+1][j]);


		 
        F_W[5][i][j+1] =((magnetic_L[i][j+1]*velocity_xL[i][j+1])+(magnetic_R[i][j+1]*velocity_xR[i][j+1]))*mu*.5 - ((bx_L[i][j+1]*velocity_yL[i][j+1])+(bx_R[i][j+1]*velocity_yR[i][j+1]))*.5-.5*(UL2[5][i][j+1]);  

        G_W[7][i+1][j] =((BX_L[i+1][j]*v_yL[i+1][j])+(BX_R[i+1][j]*v_yR[i+1][j]))*mu*.5 - ((BY_L[i+1][j]*v_xL[i+1][j])+(BY_R[i+1][j]*v_xR[i+1][j]))*.5-.5*(UR2[7][i+1][j]);  

	  
 

        F_W[6][i][j+1] = ((bz_L[i][j+1]*velocity_xL[i][j+1])+(bz_R[i][j+1]*velocity_xR[i][j+1]))*mu*.5 - ((bx_R[i][j+1]*velocity_zR[i][j+1])+(bx_L[i][j+1]*velocity_zL[i][j+1]))*.5-.5*(UL2[6][i][j+1]);

        G_W[6][i+1][j] = ((BZ_L[i+1][i]*v_yL[i+1][j])+(BZ_R[i+1][j]*v_yR[i+1][j]))*mu*.5 - ((BY_L[i+1][j]*v_zL[i+1][j])+(BY_R[i+1][j]*v_zR[i+1][j]))*.5-.5*(UR2[6][i+1][j]);
    


	  
	  
	    F_W[4][i][j+1] = ((density_L[i][j+1]*velocity_xL[i][j+1]*ethalpy_L[i][j+1]) + (density_R[i][j+1]*
                      velocity_xR[i][j+1]*ethalpy_R[i][j+1]))*.5- (.5*mu*(bx_L[i][j+1]*(bx_L[i][j+1]*velocity_xL[i][j+1]+(magnetic_L[i][j+1]*velocity_yL[i][j+1])+(bz_L[i][j+1]*velocity_zL[i][j+1]))+(bx_R[i][j+1]*(bx_R[i][j+1]*velocity_xR[i][j+1]+(magnetic_R[i][j+1]*velocity_yR[i][j+1])+(bz_R[i][j+1]*velocity_zR[i][j+1])))))-.5*(UL2[4][i][j+1]); 

  
        G_W[4][i+1][j] = ((d_L[i+1][j]*v_yL[i+1][j]*e_L[i+1][j]) + (d_R[i+1][j]*
                        v_yR[i+1][j]*e_R[i+1][j]))*.5- (.5*mu*(BY_L[i+1][j]*(BX_L[i+1][j]*v_xL[i+1][j]+(BY_L[i+1][j]*v_yL[i+1][j])+(BZ_L[i+1][j]*v_zL[i+1][j]))+(BY_R[i+1][j]*(BX_R[i+1][j]*v_xR[i+1][j]+(BY_R[i+1][j]*v_yR[i+1][j])+(BZ_R[i+1][j]*v_zR[i+1][j])))))-.5*(UR2[4][i+1][j]); 

 
        
 F_W[7][i][j+1] = -.5*(UL2[7][i][j+1]);


    G_W[5][i+1][j] = -.5*(UR2[5][i+1][j]);
        
   

  //if(i==50)
//cout<<G_W[1][i+1][j]<<"\t"<<t<<"\n";

}

}



}


void hybrid()
{


for(i=6; i<=size-8; i=i+2)
    {
	                                            
  for(j=6; j<=size-8; j=j+2)
    {

  
	  
	  for(k=0;k<8;k++)
      {
/*8
 if(fabs(density[i][j+2]-density[i][j])>1e-01)
 W1=0.2;
 else
 W1=0.0;
 
 if(fabs(density[i+2][j]-density[i][j])>1e-01)
 W2=0.2;
 else
 W2=0.0;
 
 */
 


    if(t <=0.038 )
      W1=0.22;
   else if(t>0.038 && t <=0.041)
      W1=0.22;
    else 
      W1=0.01;
   

       F1[k][i][j+1] = W1*F1[k][i][j+1]+  (1-W1)*F_W[k][i][j+1];

       G1[k][i+1][j]  = W1*G1[k][i+1][j] + (1-W1)*G_W[k][i+1][j];

             }
          }
       
       } 
   }	
	
	
void source()
{

 for(i=8; i<=size-10; i=i+2)       // need to be checked // some assumptions ? is small though
    {
   for(j=8; j<=size-10; j=j+2)
   
   {
   
   
   
 divBx[i+2][j+2]= (bx[i+2][j+2] + bx[i+2][j])/2;// + bx[i][j] + bx[i][j+2] )/4;
   
 divBy[i+2][j+2]= (magnetic[i+2][j+2] + magnetic[i][j+2])/2;// + magnetic[i][j] + magnetic[i+2][j])/4;
   
   
   
}
   }

 for(i=8; i<=size-12; i=i+2)       // need to be checked // some assumptions ? is small though
    {
   for(j=8; j<=size-12; j=j+2)
   
   {


 NetdivBx[i+2][j+2] =  (divBx[i+2][j+4]-divBx[i+2][j+2])/delta_x; 
   
 NetdivBy[i+2][j+2] =  (divBy[i+4][j+2]-divBy[i+2][j+2])/delta_x;
 

//if(i==50)
 //cout<<NetdivBy[i][j]<<"\t"<<t<<endl;

}

}


  
	 for(i=8; i<=size-12; i=i+2)       // need to be checked // some assumptions ? is small though
    {
   for(j=8; j<=size-12; j=j+2)
    {
                S[0][i+2][j+2]=0;
        
		S[1][i+2][j+2] = (NetdivBx[i+2][j+2]+ NetdivBy[i+2][j+2])*bx[i+2][j+2]; 
		
		S[2][i+2][j+2] = (NetdivBx[i+2][j+2]+ NetdivBy[i+2][j+2])*magnetic[i+2][j+2]; 
		 		
		S[3][i+2][j+2] = (NetdivBx[i+2][j+2]+ NetdivBy[i+2][j+2])*Bz[i+2][j+2]; 
		 		
		S[4][i+2][j+2] = ((bx[i+2][j+2]*velocity_x[i+2][j+2])+(magnetic[i+2][j+2]*velocity_y[i+2][j+2])+(Bz[i+2][j+2]*velocity_z[i+2][j+2]))*(NetdivBx[i+2][j+2] +NetdivBy[i+2][j+2]); 
		
		S[5][i+2][j+2] =  velocity_y[i+2][j+2]*(NetdivBx[i+2][j+2]+ NetdivBy[i+2][j+2]);

		S[6][i+2][j+2] = velocity_z[i+2][j+2]*(NetdivBx[i+2][j+2]+ NetdivBy[i+2][j+2]);
        
		S[7][i+2][j+2] = velocity_x[i+2][j+2]*(NetdivBx[i+2][j+2]+ NetdivBy[i+2][j+2]);  

  
 //  if(i==10)
//cout<<S[4][i][j]<<"\t"<<t<<"\n";
  
    }
}

return;
   }
	   
	   

	
