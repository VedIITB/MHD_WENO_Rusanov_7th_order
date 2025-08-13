#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>

#define  GAMMA                1.67
#define  max_no_variable      8
#define  L                    1.0
#define  mu                   1
#define  ghost                8
#define  cells                500
#define  faces                cells
#define  points               cells+faces+ghost+9
#define  size                 points-1

#define   gamma1     0.0286
#define   gamma2     0.343
#define   gamma3     0.5142
#define   gamma4     0.1142

#define cr1      0.0833
#define cr2      0.25
#define cr3      0.4167
#define cr4      0.5833
#define cr5      1.0833
#define cr6      1.9167
#define cr7      2.0833

#define    elipson    1e-06

 
using namespace std;


int i,j,k;
 
double  SI1,SI2,SI3,SI4,omega1,omega2,omega3,omega4,omega_1,omega_2,omega_3,omega_4,w1,w2,w3,w4,w_4,w_1,w_2,w_3,S1,S2,S3,S4,S_4,S_1,S_2,S_3;

double pressure_L[points][points],pressure_R[points][points],ethalpy_R[points][points],density_L[points][points],density_R[points][points],ethalpy_L[points][points],velocity_yR[points][points],velocity_yL[points][points];                               

double velocity_xR[points][points],velocity_xL[points][points],magnetic_L[points][points],magnetic_R[points][points],bx_R[points][points],bx_L[points][points];

double p_L[points][points],p_R[points][points],e_R[points][points],d_L[points][points],d_R[points][points],e_L[points][points],v_yR[points][points],v_yL[points][points];                               

double v_xR[points][points],v_xL[points][points],BY_L[points][points],BY_R[points][points],BX_R[points][points],BX_L[points][points], bz_R[points][points],bz_L[points][points], BZ_R[points][points], BZ_L[points][points] ;

double velocity_zL[points][points], velocity_zR[points][points] , v_zL[points][points], v_zR[points][points];


void WENO1(double u[][points]);
void WENO2(double u1[][points]);
void WENO3(double u2[][points]);
void WENO4(double u3[][points]);
void WENO5(double u4[][points]);
void WENO6(double u5[][points]); 
void WENO7(double u6[][points]); 
void WENO8(double u7[][points]);  
void WENO9(double u8[][points]);  
 
double VANLEER1(double l1 ,double l2 , double l3 );


double VANLEER1(double l1 ,double l2 , double l3 )

{

    double phi,r;

   r= (l1-l3)/((l2-l1)+1e-8);

    if(r<=0) 
     phi=0;
 
    else if(r>0)
    phi = (r+abs(r))/(1+abs(r));
 
   

    return phi;

         }



 void WENO1(double u[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
 {
  
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// for x-direction
 
 SI1 = u[i][j-6]*(547*u[i][j-6]-3882*u[i][j-4]+4642*u[i][j-2]-1854*u[i][j]) + u[i][j-4]*(7043*u[i][j-4]-17246*u[i][j-2]+7042*u[i][j])+ u[i][j-2]*(11003*u[i][j-2]-9402*u[i][j])+2107*u[i][j]*u[i][j];
 
 SI2 = u[i][j-4]*(267*u[i][j-4]-1642*u[i][j-2] + 1602*u[i][j]-494*u[i][j+2])+ u[i][j-2]*(2843*u[i][j-2]-5966*u[i][j] + 1922*u[i][j+2]) +u[i][j]*(3443*u[i][j]-2522*u[i][j+2])+ (u[i][j+2]*u[i][j+2]*547);
 
 SI3 = u[i][j-2]*(547*u[i][j-2]-2522*u[i][j] + 1922*u[i][j+2]- 494*u[i][j+4]) + u[i][j]*(3443*u[i][j]-5966*u[i][j+2]+1602*u[i][j+4]) + u[i][j+2]*(2843*u[i][j+2]-1642*u[i][j+4]) + u[i][j+4]*(267*u[i][j+4]);
 
 SI4 = u[i][j]*(2107*u[i][j]-9402*u[i][j+2] + 7042*u[i][j+4]- 1854*u[i][j+6]) + u[i][j+2]*(11003*u[i][j+2]-17246*u[i][j+4]+4642*u[i][j+6]) + u[i][j+4]*(7043*u[i][j+4]-3882*u[i][j+6]) + u[i][j+6]*(547*u[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u[i][j-6] + (cr5)*u[i][j-4] - (cr6)*u[i][j-2]+ (cr7)*u[i][j]; 
 	
 S2 = (cr1)*u[i][j-4] - (cr3)*u[i][j-2] + (cr5)*u[i][j] + (cr2)*u[i][j+2]; 
 
 S3 =  (-cr1)*u[i][j-2] + (cr4)*u[i][j] + (cr4)*u[i][j+2]- (cr1)*u[i][j+4]; 
 
 S4 =  (-cr3)*u[i][j+4] + (cr2)*u[i][j] + (cr5)*u[i][j+2]+ (cr1)*u[i][j+6]; 
 	
 	
 S_4 = (-cr2)*u[i][j+6] + (cr5)*u[i][j+4] - (cr6)*u[i][j+2]+ (cr7)*u[i][j];
 
 S_3 =  cr1*u[i][j+4] - cr3*u[i][j+2] + cr5*u[i][j] + cr2*u[i][j-2]; 
 
 S_2 = -cr1*u[i][j+2] + cr4*u[i][j] + cr4*u[i][j-2] - cr1*u[i][j-4] ; 
 
 S_1 =  cr2*u[i][j]+ cr5*u[i][j-2] - cr3*u[i][j-4]+ cr1*u[i][j-6];
 
  
  density_L[i][j+1]=(S1*w1)+(S2*w2)+(S3*w3)+ S4*w4;	
 	
  density_R[i][j-1]=(S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
 	
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction reconstruction
  
  
 SI1 = u[i-6][j]*(547*u[i-6][j]-3882*u[i-4][j]+4642*u[i-2][j]-1854*u[i][j]) + u[i-4][j]*(7043*u[i-4][j]-17246*u[i-2][j]+7042*u[i][j])+ u[i-2][j]*(11003*u[i-2][j]-9402*u[i][j])+2107*u[i][j]*u[i][j];
 
 SI2 = u[i-4][j]*(267*u[i-4][j]-1642*u[i-2][j] + 1602*u[i][j]-494*u[i+2][j])+ u[i-2][j]*(2843*u[i-2][j]-5966*u[i][j] + 1922*u[i+2][j]) +u[i][j]*(3443*u[i][j]-2522*u[i+2][j])+ (u[i+2][j]*u[i+2][j]*547);
 
 SI3 = u[i-2][j]*(547*u[i-2][j]-2522*u[i][j] + 1922*u[i+2][j]- 494*u[i+4][j]) + u[i][j]*(3443*u[i][j]-5966*u[i+2][j]+1602*u[i+4][j]) + u[i+2][j]*(2843*u[i+2][j]-1642*u[i+4][j]) + u[i+4][j]*(267*u[i+4][j]);
 
 SI4 = u[i][j]*(2107*u[i][j]-9402*u[i+2][j] + 7042*u[i+4][j]- 1854*u[i+6][j]) + u[i+2][j]*(11003*u[i+2][j]-17246*u[i+4][j]+4642*u[i+6][j]) + u[i+4][j]*(7043*u[i+4][j]-3882*u[i+6][j]) + u[i+6][j]*(547*u[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u[i-6][j] + (cr5)*u[i-4][j] - (cr6)*u[i-2][j]+ (cr7)*u[i][j]; 
 	
 S2 = (cr1)*u[i-4][j] - (cr3)*u[i-2][j] + (cr5)*u[i][j] + (cr2)*u[i+2][j]; 
 
 S3 =  (-cr1)*u[i-2][j] + cr4*u[i][j] + cr4*u[i+2][j]- cr1*u[i+4][j]; 
 
 S4 =  (-cr3)*u[i+4][j] + cr2*u[i][j] + cr5*u[i+2][j]+ cr1*u[i+6][j]; 
 	
 	
 S_4 = -cr2*u[i+6][j] + cr5*u[i+4][j] - cr6*u[i+2][j]+ cr7*u[i][j];
 
 S_3 =  cr1*u[i+4][j] - cr3*u[i+2][j] + cr5*u[i][j] + cr2*u[i-2][j]; 
 
 S_2 = -cr1*u[i+2][j] + cr4*u[i][j] + cr4*u[i-2][j] - cr1*u[i-4][j] ; 
 
 S_1 =  cr2*u[i][j]+ cr5*u[i-2][j] - cr3*u[i-4][j]+ cr1*u[i-6][j];
 
  
  
  d_L[i+1][j]=(S1*w1)+(S2*w2)+(S3*w3)+ S4*w4;	
 	
  d_R[i-1][j]=(S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
  
 	 
	  //if(i==50)
	
 	 }
}

return;
}



 void WENO2(double u1[][points])
 {


for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
  SI1 = u1[i][j-6]*(547*u1[i][j-6]-3882*u1[i][j-4]+4642*u1[i][j-2]-1854*u1[i][j]) + u1[i][j-4]*(7043*u1[i][j-4]-17246*u1[i][j-2]+7042*u1[i][j])+ u1[i][j-2]*(11003*u1[i][j-2]-9402*u1[i][j])+2107*u1[i][j]*u1[i][j];
 
 SI2 = u1[i][j-4]*(267*u1[i][j-4]-1642*u1[i][j-2] + 1602*u1[i][j]-494*u1[i][j+2])+ u1[i][j-2]*(2843*u1[i][j-2]-5966*u1[i][j] + 1922*u1[i][j+2]) +u1[i][j]*(3443*u1[i][j]-2522*u1[i][j+2])+ (u1[i][j+2]*u1[i][j+2]*547);
 
 SI3 = u1[i][j-2]*(547*u1[i][j-2]-2522*u1[i][j] + 1922*u1[i][j+2]- 494*u1[i][j+4]) + u1[i][j]*(3443*u1[i][j]-5966*u1[i][j+2]+1602*u1[i][j+4]) + u1[i][j+2]*(2843*u1[i][j+2]-1642*u1[i][j+4]) + u1[i][j+4]*(267*u1[i][j+4]);
 
 SI4 = u1[i][j]*(2107*u1[i][j]-9402*u1[i][j+2] + 7042*u1[i][j+4]- 1854*u1[i][j+6]) + u1[i][j+2]*(11003*u1[i][j+2]-17246*u1[i][j+4]+4642*u1[i][j+6]) + u1[i][j+4]*(7043*u1[i][j+4]-3882*u1[i][j+6]) + u1[i][j+6]*(547*u1[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u1[i][j-6] + (cr5)*u1[i][j-4] - (cr6)*u1[i][j-2]+ (cr7)*u1[i][j]; 
 	
 S2 = (cr1)*u1[i][j-4] - (cr3)*u1[i][j-2] + (cr5)*u1[i][j] + (cr2)*u1[i][j+2]; 
 
 S3 =  (-cr1)*u1[i][j-2] + cr4*u1[i][j] + cr4*u1[i][j+2]- cr1*u1[i][j+4]; 
 
 S4 =  (-cr3)*u1[i][j+4] + cr2*u1[i][j] + cr5*u1[i][j+2]+ cr1*u1[i][j+6]; 
 	
 	
 S_4 = -cr2*u1[i][j+6] + cr5*u1[i][j+4] - cr6*u1[i][j+2]+ cr7*u1[i][j];
 
 S_3 =  cr1*u1[i][j+4] - cr3*u1[i][j+2] + cr5*u1[i][j] + cr2*u1[i][j-2]; 
 
 S_2 = -cr1*u1[i][j+2] + cr4*u1[i][j] + cr4*u1[i][j-2] - cr1*u1[i][j-4] ; 
 
 S_1 =  cr2*u1[i][j]+ cr5*u1[i][j-2] - cr3*u1[i][j-4]+ cr1*u1[i][j-6];
 

 
  pressure_L[i][j+1]=S1*w1+S2*w2+S3*w3 +S4*w4;	
 	
  pressure_R[i][j-1]=S_1*w_1+S_2*w_2+S_3*w_3+S_4*w_4;	
 	
 	
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
 SI1 = u1[i-6][j]*(547*u1[i-6][j]-3882*u1[i-4][j]+4642*u1[i-2][j]-1854*u1[i][j]) + u1[i-4][j]*(7043*u1[i-4][j]-17246*u1[i-2][j]+7042*u1[i][j])+ u1[i-2][j]*(11003*u1[i-2][j]-9402*u1[i][j])+2107*u1[i][j]*u1[i][j];
 
 SI2 = u1[i-4][j]*(267*u1[i-4][j]-1642*u1[i-2][j] + 1602*u1[i][j]-494*u1[i+2][j])+ u1[i-2][j]*(2843*u1[i-2][j]-5966*u1[i][j] + 1922*u1[i+2][j]) +u1[i][j]*(3443*u1[i][j]-2522*u1[i+2][j])+ (u1[i+2][j]*u1[i+2][j]*547);
 
 SI3 = u1[i-2][j]*(547*u1[i-2][j]-2522*u1[i][j] + 1922*u1[i+2][j]- 494*u1[i+4][j]) + u1[i][j]*(3443*u1[i][j]-5966*u1[i+2][j]+1602*u1[i+4][j]) + u1[i+2][j]*(2843*u1[i+2][j]-1642*u1[i+4][j]) + u1[i+4][j]*(267*u1[i+4][j]);
 
 SI4 = u1[i][j]*(2107*u1[i][j]-9402*u1[i+2][j] + 7042*u1[i+4][j]- 1854*u1[i+6][j]) + u1[i+2][j]*(11003*u1[i+2][j]-17246*u1[i+4][j]+4642*u1[i+6][j]) + u1[i+4][j]*(7043*u1[i+4][j]-3882*u1[i+6][j]) + u1[i+6][j]*(547*u1[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u1[i-6][j] + (cr5)*u1[i-4][j] - (cr6)*u1[i-2][j]+ (cr7)*u1[i][j]; 
 	
 S2 = (cr1)*u1[i-4][j] - (cr3)*u1[i-2][j] + (cr5)*u1[i][j] + (cr2)*u1[i+2][j]; 
 
 S3 =  (-cr1)*u1[i-2][j] + cr4*u1[i][j] + cr4*u1[i+2][j]- cr1*u1[i+4][j]; 
 
 S4 =  (-cr3)*u1[i+4][j] + cr2*u1[i][j] + cr5*u1[i+2][j]+ cr1*u1[i+6][j]; 
 	
 	
 S_4 = -cr2*u1[i+6][j] + cr5*u1[i+4][j] - cr6*u1[i+2][j]+ cr7*u1[i][j];
 
 S_3 =  cr1*u1[i+4][j] - cr3*u1[i+2][j] + cr5*u1[i][j] + cr2*u1[i-2][j]; 
 
 S_2 = -cr1*u1[i+2][j] + cr4*u1[i][j] + cr4*u1[i-2][j] - cr1*u1[i-4][j] ; 
 
 S_1 =  cr2*u1[i][j]+ cr5*u1[i-2][j] - cr3*u1[i-4][j]+ cr1*u1[i-6][j];
 
 
  
 	
 p_L[i+1][j]=S1*w1+S2*w2+S3*w3 + S4*w4;	
 	
  p_R[i-1][j]=S_1*w_1+S_2*w_2+S_3*w_3 +S_4*w_4;	
 	
 	
 	  //cout<< p_L[i+1][j]<<"\t"<<i<<endl;
 	
	  } 
  
  //return;
  }	 
}

void WENO3(double u2[][points])
{

for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
{
    
 // if(i==10)
  //cin>>j;
 
 SI1 = u2[i][j-6]*(547*u2[i][j-6]-3882*u2[i][j-4]+4642*u2[i][j-2]-1854*u2[i][j]) + u2[i][j-4]*(7043*u2[i][j-4]-17246*u2[i][j-2]+7042*u2[i][j])+ u2[i][j-2]*(11003*u2[i][j-2]-9402*u2[i][j])+2107*u2[i][j]*u2[i][j];
 
 SI2 = u2[i][j-4]*(267*u2[i][j-4]-1642*u2[i][j-2] + 1602*u2[i][j]-494*u2[i][j+2])+ u2[i][j-2]*(2843*u2[i][j-2]-5966*u2[i][j] + 1922*u2[i][j+2]) +u2[i][j]*(3443*u2[i][j]-2522*u2[i][j+2])+ (u2[i][j+2]*u2[i][j+2]*547);
 
 SI3 = u2[i][j-2]*(547*u2[i][j-2]-2522*u2[i][j] + 1922*u2[i][j+2]- 494*u2[i][j+4]) + u2[i][j]*(3443*u2[i][j]-5966*u2[i][j+2]+1602*u2[i][j+4]) + u2[i][j+2]*(2843*u2[i][j+2]-1642*u2[i][j+4]) + u2[i][j+4]*(267*u2[i][j+4]);
 
 SI4 = u2[i][j]*(2107*u2[i][j]-9402*u2[i][j+2] + 7042*u2[i][j+4]- 1854*u2[i][j+6]) + u2[i][j+2]*(11003*u2[i][j+2]-17246*u2[i][j+4]+4642*u2[i][j+6]) + u2[i][j+4]*(7043*u2[i][j+4]-3882*u2[i][j+6]) + u2[i][j+6]*(547*u2[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u2[i][j-6] + (cr5)*u2[i][j-4] - (cr6)*u2[i][j-2]+ (cr7)*u2[i][j]; 
 	
 S2 = (cr1)*u2[i][j-4] - (cr3)*u2[i][j-2] + (cr5)*u2[i][j] + (cr2)*u2[i][j+2]; 
 
 S3 =  (-cr1)*u2[i][j-2] + cr4*u2[i][j] + cr4*u2[i][j+2]- cr1*u2[i][j+4]; 
 
 S4 =  (-cr3)*u2[i][j+4] + cr2*u2[i][j] + cr5*u2[i][j+2]+ cr1*u2[i][j+6]; 
 	
 	
 S_4 = -cr2*u2[i][j+6] + cr5*u2[i][j+4] - cr6*u2[i][j+2]+ cr7*u2[i][j];
 
 S_3 =  cr1*u2[i][j+4] - cr3*u2[i][j+2] + cr5*u2[i][j] + cr2*u2[i][j-2]; 
 
 S_2 = -cr1*u2[i][j+2] + cr4*u2[i][j] + cr4*u2[i][j-2] - cr1*u2[i][j-4] ; 
 
 S_1 =  cr2*u2[i][j]+ cr5*u2[i][j-2] - cr3*u2[i][j-4]+ cr1*u2[i][j-6];
 
                                                          
  velocity_yL[i][j+1] =S1*w1+S2*w2+S3*w3 + S4*w4;	
 	
  velocity_yR[i][j-1] =S_1*w_1+S_2*w_2+S_3*w_3+ S_4*w_4;	
 	


 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
 SI1 = u2[i-6][j]*(547*u2[i-6][j]-3882*u2[i-4][j]+4642*u2[i-2][j]-1854*u2[i][j]) + u2[i-4][j]*(7043*u2[i-4][j]-17246*u2[i-2][j]+7042*u2[i][j])+ u2[i-2][j]*(11003*u2[i-2][j]-9402*u2[i][j])+2107*u2[i][j]*u2[i][j];
 
 SI2 = u2[i-4][j]*(267*u2[i-4][j]-1642*u2[i-2][j] + 1602*u2[i][j]-494*u2[i+2][j])+ u2[i-2][j]*(2843*u2[i-2][j]-5966*u2[i][j] + 1922*u2[i+2][j]) +u2[i][j]*(3443*u2[i][j]-2522*u2[i+2][j])+ (u2[i+2][j]*u2[i+2][j]*547);
 
 SI3 = u2[i-2][j]*(547*u2[i-2][j]-2522*u2[i][j] + 1922*u2[i+2][j]- 494*u2[i+4][j]) + u2[i][j]*(3443*u2[i][j]-5966*u2[i+2][j]+1602*u2[i+4][j]) + u2[i+2][j]*(2843*u2[i+2][j]-1642*u2[i+4][j]) + u2[i+4][j]*(267*u2[i+4][j]);
 
 SI4 = u2[i][j]*(2107*u2[i][j]-9402*u2[i+2][j] + 7042*u2[i+4][j]- 1854*u2[i+6][j]) + u2[i+2][j]*(11003*u2[i+2][j]-17246*u2[i+4][j]+4642*u2[i+6][j]) + u2[i+4][j]*(7043*u2[i+4][j]-3882*u2[i+6][j]) + u2[i+6][j]*(547*u2[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u2[i-6][j] + (cr5)*u2[i-4][j] - (cr6)*u2[i-2][j]+ (cr7)*u2[i][j]; 
 	
 S2 = (cr1)*u2[i-4][j] - (cr3)*u2[i-2][j] + (cr5)*u2[i][j] + (cr2)*u2[i+2][j]; 
 
 S3 =  (-cr1)*u2[i-2][j] + cr4*u2[i][j] + cr4*u2[i+2][j]- cr1*u2[i+4][j]; 
 
 S4 =  (-cr3)*u2[i+4][j] + cr2*u2[i][j] + cr5*u2[i+2][j]+ cr1*u2[i+6][j]; 
 	
 	
 S_4 = -cr2*u2[i+6][j] + cr5*u2[i+4][j] - cr6*u2[i+2][j]+ cr7*u2[i][j];
 
 S_3 =  cr1*u2[i+4][j] - cr3*u2[i+2][j] + cr5*u2[i][j] + cr2*u2[i-2][j]; 
 
 S_2 = -cr1*u2[i+2][j] + cr4*u2[i][j] + cr4*u2[i-2][j] - cr1*u2[i-4][j] ; 
 
 S_1 =  cr2*u2[i][j]+ cr5*u2[i-2][j] - cr3*u2[i-4][j]+ cr1*u2[i-6][j];
  
  


v_yL[i+1][j] =S1*w1+S2*w2+S3*w3+S4*w4;	
 	
  v_yR[i-1][j] =S_1*w_1+S_2*w_2+S_3*w_3 + S_4*w_4;	
 
 
 	 }
 	// return;
}
}


void WENO4(double u3[][points])
{

for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
{
    
 // if(i==10)
  //cin>>j;
 
 
 SI1 = u3[i][j-6]*(547*u3[i][j-6]-3882*u3[i][j-4]+4642*u3[i][j-2]-1854*u3[i][j]) + u3[i][j-4]*(7043*u3[i][j-4]-17246*u3[i][j-2]+7042*u3[i][j])+ u3[i][j-2]*(11003*u3[i][j-2]-9402*u3[i][j])+2107*u3[i][j]*u3[i][j];
 
 SI2 = u3[i][j-4]*(267*u3[i][j-4]-1642*u3[i][j-2] + 1602*u3[i][j]-494*u3[i][j+2])+ u3[i][j-2]*(2843*u3[i][j-2]-5966*u3[i][j] + 1922*u3[i][j+2]) +u3[i][j]*(3443*u3[i][j]-2522*u3[i][j+2])+ (u3[i][j+2]*u3[i][j+2]*547);
 
 SI3 = u3[i][j-2]*(547*u3[i][j-2]-2522*u3[i][j] + 1922*u3[i][j+2]- 494*u3[i][j+4]) + u3[i][j]*(3443*u3[i][j]-5966*u3[i][j+2]+1602*u3[i][j+4]) + u3[i][j+2]*(2843*u3[i][j+2]-1642*u3[i][j+4]) + u3[i][j+4]*(267*u3[i][j+4]);
 
 SI4 = u3[i][j]*(2107*u3[i][j]-9402*u3[i][j+2] + 7042*u3[i][j+4]- 1854*u3[i][j+6]) + u3[i][j+2]*(11003*u3[i][j+2]-17246*u3[i][j+4]+4642*u3[i][j+6]) + u3[i][j+4]*(7043*u3[i][j+4]-3882*u3[i][j+6]) + u3[i][j+6]*(547*u3[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u3[i][j-6] + (cr5)*u3[i][j-4] - (cr6)*u3[i][j-2]+ (cr7)*u3[i][j]; 
 	
 S2 = (cr1)*u3[i][j-4] - (cr3)*u3[i][j-2] + (cr5)*u3[i][j] + (cr2)*u3[i][j+2]; 
 
 S3 =  (-cr1)*u3[i][j-2] + cr4*u3[i][j] + cr4*u3[i][j+2]- cr1*u3[i][j+4]; 
 
 S4 =  (-cr3)*u3[i][j+4] + cr2*u3[i][j] + cr5*u3[i][j+2]+ cr1*u3[i][j+6]; 
 	
 	
 S_4 = -cr2*u3[i][j+6] + cr5*u3[i][j+4] - cr6*u3[i][j+2]+ cr7*u3[i][j];
 
 S_3 =  cr1*u3[i][j+4] - cr3*u3[i][j+2] + cr5*u3[i][j] + cr2*u3[i][j-2]; 
 
 S_2 = -cr1*u3[i][j+2] + cr4*u3[i][j] + cr4*u3[i][j-2] - cr1*u3[i][j-4] ; 
 
 S_1 =  cr2*u3[i][j]+ cr5*u3[i][j-2] - cr3*u3[i][j-4]+ cr1*u3[i][j-6];
 
 

                                                           
  ethalpy_L[i][j+1] = S1*w1+S2*w2+S3*w3 +S4*w4;	
 	
  ethalpy_R[i][j-1] = S_1*w_1+S_2*w_2+S_3*w_3 + S_4*w_4;	
 	
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
 SI1 = u3[i-6][j]*(547*u3[i-6][j]-3882*u3[i-4][j]+4642*u3[i-2][j]-1854*u3[i][j]) + u3[i-4][j]*(7043*u3[i-4][j]-17246*u3[i-2][j]+7042*u3[i][j])+ u3[i-2][j]*(11003*u3[i-2][j]-9402*u3[i][j])+2107*u3[i][j]*u3[i][j];
 
 SI2 = u3[i-4][j]*(267*u3[i-4][j]-1642*u3[i-2][j] + 1602*u3[i][j]-494*u3[i+2][j])+ u3[i-2][j]*(2843*u3[i-2][j]-5966*u3[i][j] + 1922*u3[i+2][j]) +u3[i][j]*(3443*u3[i][j]-2522*u3[i+2][j])+ (u3[i+2][j]*u3[i+2][j]*547);
 
 SI3 = u3[i-2][j]*(547*u3[i-2][j]-2522*u3[i][j] + 1922*u3[i+2][j]- 494*u3[i+4][j]) + u3[i][j]*(3443*u3[i][j]-5966*u3[i+2][j]+1602*u3[i+4][j]) + u3[i+2][j]*(2843*u3[i+2][j]-1642*u3[i+4][j]) + u3[i+4][j]*(267*u3[i+4][j]);
 
 SI4 = u3[i][j]*(2107*u3[i][j]-9402*u3[i+2][j] + 7042*u3[i+4][j]- 1854*u3[i+6][j]) + u3[i+2][j]*(11003*u3[i+2][j]-17246*u3[i+4][j]+4642*u3[i+6][j]) + u3[i+4][j]*(7043*u3[i+4][j]-3882*u3[i+6][j]) + u3[i+6][j]*(547*u3[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u3[i-6][j] + (cr5)*u3[i-4][j] - (cr6)*u3[i-2][j]+ (cr7)*u3[i][j]; 
 	
 S2 = (cr1)*u3[i-4][j] - (cr3)*u3[i-2][j] + (cr5)*u3[i][j] + (cr2)*u3[i+2][j]; 
 
 S3 =  (-cr1)*u3[i-2][j] + cr4*u3[i][j] + cr4*u3[i+2][j]- cr1*u3[i+4][j]; 
 
 S4 =  (-cr3)*u3[i+4][j] + cr2*u3[i][j] + cr5*u3[i+2][j]+ cr1*u3[i+6][j]; 
 	
 	
 S_4 = -cr2*u3[i+6][j] + cr5*u3[i+4][j] - cr6*u3[i+2][j]+ cr7*u3[i][j];
 
 S_3 =  cr1*u3[i+4][j] - cr3*u3[i+2][j] + cr5*u3[i][j] + cr2*u3[i-2][j]; 
 
 S_2 = -cr1*u3[i+2][j] + cr4*u3[i][j] + cr4*u3[i-2][j] - cr1*u3[i-4][j] ; 
 
 S_1 =  cr2*u3[i][j]+ cr5*u3[i-2][j] - cr3*u3[i-4][j]+ cr1*u3[i-6][j];

	
e_L[i+1][j] = S1*w1+S2*w2+S3*w3+S4*w4;	
 	
  e_R[i-1][j] = S_1*w_1+S_2*w_2+S_3*w_3 + S_4*w_4;	
 
 	 
 	 }
 //return;
}
}
void WENO5(double u4[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
SI1 = u4[i][j-6]*(547*u4[i][j-6]-3882*u4[i][j-4]+4642*u4[i][j-2]-1854*u4[i][j]) + u4[i][j-4]*(7043*u4[i][j-4]-17246*u4[i][j-2]+7042*u4[i][j])+ u4[i][j-2]*(11003*u4[i][j-2]-9402*u4[i][j])+2107*u4[i][j]*u4[i][j];
 
 SI2 = u4[i][j-4]*(267*u4[i][j-4]-1642*u4[i][j-2] + 1602*u4[i][j]-494*u4[i][j+2])+ u4[i][j-2]*(2843*u4[i][j-2]-5966*u4[i][j] + 1922*u4[i][j+2]) +u4[i][j]*(3443*u4[i][j]-2522*u4[i][j+2])+ (u4[i][j+2]*u4[i][j+2]*547);
 
 SI3 = u4[i][j-2]*(547*u4[i][j-2]-2522*u4[i][j] + 1922*u4[i][j+2]- 494*u4[i][j+4]) + u4[i][j]*(3443*u4[i][j]-5966*u4[i][j+2]+1602*u4[i][j+4]) + u4[i][j+2]*(2843*u4[i][j+2]-1642*u4[i][j+4]) + u4[i][j+4]*(267*u4[i][j+4]);
 
 SI4 = u4[i][j]*(2107*u4[i][j]-9402*u4[i][j+2] + 7042*u4[i][j+4]- 1854*u4[i][j+6]) + u4[i][j+2]*(11003*u4[i][j+2]-17246*u4[i][j+4]+4642*u4[i][j+6]) + u4[i][j+4]*(7043*u4[i][j+4]-3882*u4[i][j+6]) + u4[i][j+6]*(547*u4[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u4[i][j-6] + (cr5)*u4[i][j-4] - (cr6)*u4[i][j-2]+ (cr7)*u4[i][j]; 
 	
 S2 = (cr1)*u4[i][j-4] - (cr3)*u4[i][j-2] + (cr5)*u4[i][j] + (cr2)*u4[i][j+2]; 
 
 S3 =  (-cr1)*u4[i][j-2] + cr4*u4[i][j] + cr4*u4[i][j+2]- cr1*u4[i][j+4]; 
 
 S4 =  (-cr3)*u4[i][j+4] + cr2*u4[i][j] + cr5*u4[i][j+2]+ cr1*u4[i][j+6]; 
 	
 	
 S_4 = -cr2*u4[i][j+6] + cr5*u4[i][j+4] - cr6*u4[i][j+2]+ cr7*u4[i][j];
 
 S_3 =  cr1*u4[i][j+4] - cr3*u4[i][j+2] + cr5*u4[i][j] + cr2*u4[i][j-2]; 
 
 S_2 = -cr1*u4[i][j+2] + cr4*u4[i][j] + cr4*u4[i][j-2] - cr1*u4[i][j-4] ; 
 
 S_1 =  cr2*u4[i][j]+ cr5*u4[i][j-2] - cr3*u4[i][j-4]+ cr1*u4[i][j-6];
 
   
 
  velocity_xL[i][j+1] = (S1*w1)+(S2*w2)+(S3*w3)+S4*w4;	
 	
  velocity_xR[i][j-1] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
 	
    
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = u4[i-6][j]*(547*u4[i-6][j]-3882*u4[i-4][j]+4642*u4[i-2][j]-1854*u4[i][j]) + u4[i-4][j]*(7043*u4[i-4][j]-17246*u4[i-2][j]+7042*u4[i][j])+ u4[i-2][j]*(11003*u4[i-2][j]-9402*u4[i][j])+2107*u4[i][j]*u4[i][j];
 
 SI2 = u4[i-4][j]*(267*u4[i-4][j]-1642*u4[i-2][j] + 1602*u4[i][j]-494*u4[i+2][j])+ u4[i-2][j]*(2843*u4[i-2][j]-5966*u4[i][j] + 1922*u4[i+2][j]) +u4[i][j]*(3443*u4[i][j]-2522*u4[i+2][j])+ (u4[i+2][j]*u4[i+2][j]*547);
 
 SI3 = u4[i-2][j]*(547*u4[i-2][j]-2522*u4[i][j] + 1922*u4[i+2][j]- 494*u4[i+4][j]) + u4[i][j]*(3443*u4[i][j]-5966*u4[i+2][j]+1602*u4[i+4][j]) + u4[i+2][j]*(2843*u4[i+2][j]-1642*u4[i+4][j]) + u4[i+4][j]*(267*u4[i+4][j]);
 
 SI4 = u4[i][j]*(2107*u4[i][j]-9402*u4[i+2][j] + 7042*u4[i+4][j]- 1854*u4[i+6][j]) + u4[i+2][j]*(11003*u4[i+2][j]-17246*u4[i+4][j]+4642*u4[i+6][j]) + u4[i+4][j]*(7043*u4[i+4][j]-3882*u4[i+6][j]) + u4[i+6][j]*(547*u4[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u4[i-6][j] + (cr5)*u4[i-4][j] - (cr6)*u4[i-2][j]+ (cr7)*u4[i][j]; 
 	
 S2 = (cr1)*u4[i-4][j] - (cr3)*u4[i-2][j] + (cr5)*u4[i][j] + (cr2)*u4[i+2][j]; 
 
 S3 =  (-cr1)*u4[i-2][j] + cr4*u4[i][j] + cr4*u4[i+2][j]- cr1*u4[i+4][j]; 
 
 S4 =  (-cr3)*u4[i+4][j] + cr2*u4[i][j] + cr5*u4[i+2][j]+ cr1*u4[i+6][j]; 
 	
 	
 S_4 = -cr2*u4[i+6][j] + cr5*u4[i+4][j] - cr6*u4[i+2][j]+ cr7*u4[i][j];
 
 S_3 =  cr1*u4[i+4][j] - cr3*u4[i+2][j] + cr5*u4[i][j] + cr2*u4[i-2][j]; 
 
 S_2 = -cr1*u4[i+2][j] + cr4*u4[i][j] + cr4*u4[i-2][j] - cr1*u4[i-4][j] ; 
 
 S_1 =  cr2*u4[i][j]+ cr5*u4[i-2][j] - cr3*u4[i-4][j]+ cr1*u4[i-6][j];


 
 v_xL[i+1][j] = (S1*w1)+(S2*w2)+(S3*w3)+S4*w4;	
 	
  v_xR[i-1][j] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
 	
 	 
 	 }
//return;

}
}
 
 void WENO6(double u5[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
SI1 = u5[i][j-6]*(547*u5[i][j-6]-3882*u5[i][j-4]+4642*u5[i][j-2]-1854*u5[i][j]) + u5[i][j-4]*(7043*u5[i][j-4]-17246*u5[i][j-2]+7042*u5[i][j])+ u5[i][j-2]*(11003*u5[i][j-2]-9402*u5[i][j])+2107*u5[i][j]*u5[i][j];
 
 SI2 = u5[i][j-4]*(267*u5[i][j-4]-1642*u5[i][j-2] + 1602*u5[i][j]-494*u5[i][j+2])+ u5[i][j-2]*(2843*u5[i][j-2]-5966*u5[i][j] + 1922*u5[i][j+2]) +u5[i][j]*(3443*u5[i][j]-2522*u5[i][j+2])+ (u5[i][j+2]*u5[i][j+2]*547);
 
 SI3 = u5[i][j-2]*(547*u5[i][j-2]-2522*u5[i][j] + 1922*u5[i][j+2]- 494*u5[i][j+4]) + u5[i][j]*(3443*u5[i][j]-5966*u5[i][j+2]+1602*u5[i][j+4]) + u5[i][j+2]*(2843*u5[i][j+2]-1642*u5[i][j+4]) + u5[i][j+4]*(267*u5[i][j+4]);
 
 SI4 = u5[i][j]*(2107*u5[i][j]-9402*u5[i][j+2] + 7042*u5[i][j+4]- 1854*u5[i][j+6]) + u5[i][j+2]*(11003*u5[i][j+2]-17246*u5[i][j+4]+4642*u5[i][j+6]) + u5[i][j+4]*(7043*u5[i][j+4]-3882*u5[i][j+6]) + u5[i][j+6]*(547*u5[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u5[i][j-6] + (cr5)*u5[i][j-4] - (cr6)*u5[i][j-2]+ (cr7)*u5[i][j]; 
 	
 S2 = (cr1)*u5[i][j-4] - (cr3)*u5[i][j-2] + (cr5)*u5[i][j] + (cr2)*u5[i][j+2]; 
 
 S3 =  (-cr1)*u5[i][j-2] + cr4*u5[i][j] + cr4*u5[i][j+2]- cr1*u5[i][j+4]; 
 
 S4 =  (-cr3)*u5[i][j+4] + cr2*u5[i][j] + cr5*u5[i][j+2]+ cr1*u5[i][j+6]; 
 	
 	
 S_4 = -cr2*u5[i][j+6] + cr5*u5[i][j+4] - cr6*u5[i][j+2]+ cr7*u5[i][j];
 
 S_3 =  cr1*u5[i][j+4] - cr3*u5[i][j+2] + cr5*u5[i][j] + cr2*u5[i][j-2]; 
 
 S_2 = -cr1*u5[i][j+2] + cr4*u5[i][j] + cr4*u5[i][j-2] - cr1*u5[i][j-4] ; 
 
 S_1 =  cr2*u5[i][j]+ cr5*u5[i][j-2] - cr3*u5[i][j-4]+ cr1*u5[i][j-6];
 
 
  magnetic_L[i][j+1] = (S1*w1)+(S2*w2)+(S3*w3)+S4*w4;	
 	
  magnetic_R[i][j-1] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
 	
    
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = u5[i-6][j]*(547*u5[i-6][j]-3882*u5[i-4][j]+4642*u5[i-2][j]-1854*u5[i][j]) + u5[i-4][j]*(7043*u5[i-4][j]-17246*u5[i-2][j]+7042*u5[i][j])+ u5[i-2][j]*(11003*u5[i-2][j]-9402*u5[i][j])+2107*u5[i][j]*u5[i][j];
 
 SI2 = u5[i-4][j]*(267*u5[i-4][j]-1642*u5[i-2][j] + 1602*u5[i][j]-494*u5[i+2][j])+ u5[i-2][j]*(2843*u5[i-2][j]-5966*u5[i][j] + 1922*u5[i+2][j]) +u5[i][j]*(3443*u5[i][j]-2522*u5[i+2][j])+ (u5[i+2][j]*u5[i+2][j]*547);
 
 SI3 = u5[i-2][j]*(547*u5[i-2][j]-2522*u5[i][j] + 1922*u5[i+2][j]- 494*u5[i+4][j]) + u5[i][j]*(3443*u5[i][j]-5966*u5[i+2][j]+1602*u5[i+4][j]) + u5[i+2][j]*(2843*u5[i+2][j]-1642*u5[i+4][j]) + u5[i+4][j]*(267*u5[i+4][j]);
 
 SI4 = u5[i][j]*(2107*u5[i][j]-9402*u5[i+2][j] + 7042*u5[i+4][j]- 1854*u5[i+6][j]) + u5[i+2][j]*(11003*u5[i+2][j]-17246*u5[i+4][j]+4642*u5[i+6][j]) + u5[i+4][j]*(7043*u5[i+4][j]-3882*u5[i+6][j]) + u5[i+6][j]*(547*u5[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u5[i-6][j] + (cr5)*u5[i-4][j] - (cr6)*u5[i-2][j]+ (cr7)*u5[i][j]; 
 	
 S2 = (cr1)*u5[i-4][j] - (cr3)*u5[i-2][j] + (cr5)*u5[i][j] + (cr2)*u5[i+2][j]; 
 
 S3 =  (-cr1)*u5[i-2][j] + cr4*u5[i][j] + cr4*u5[i+2][j]- cr1*u5[i+4][j]; 
 
 S4 =  (-cr3)*u5[i+4][j] + cr2*u5[i][j] + cr5*u5[i+2][j]+ cr1*u5[i+6][j]; 
 	
 	
 S_4 = -cr2*u5[i+6][j] + cr5*u5[i+4][j] - cr6*u5[i+2][j]+ cr7*u5[i][j];
 
 S_3 =  cr1*u5[i+4][j] - cr3*u5[i+2][j] + cr5*u5[i][j] + cr2*u5[i-2][j]; 
 
 S_2 = -cr1*u5[i+2][j] + cr4*u5[i][j] + cr4*u5[i-2][j] - cr1*u5[i-4][j] ; 
 
 S_1 =  cr2*u5[i][j]+ cr5*u5[i-2][j] - cr3*u5[i-4][j]+ cr1*u5[i-6][j];


  BY_L[i+1][j] = (S1*w1)+(S2*w2)+(S3*w3)+S4*w4;	
 	
 BY_R[i-1][j] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3) + S_4*w_4;	
 	
 	 
 	 }


}
}
void WENO7(double u6[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
SI1 = u6[i][j-6]*(547*u6[i][j-6]-3882*u6[i][j-4]+4642*u6[i][j-2]-1854*u6[i][j]) + u6[i][j-4]*(7043*u6[i][j-4]-17246*u6[i][j-2]+7042*u6[i][j])+ u6[i][j-2]*(11003*u6[i][j-2]-9402*u6[i][j])+2107*u6[i][j]*u6[i][j];
 
 SI2 = u6[i][j-4]*(267*u6[i][j-4]-1642*u6[i][j-2] + 1602*u6[i][j]-494*u6[i][j+2])+ u6[i][j-2]*(2843*u6[i][j-2]-5966*u6[i][j] + 1922*u6[i][j+2]) +u6[i][j]*(3443*u6[i][j]-2522*u6[i][j+2])+ (u6[i][j+2]*u6[i][j+2]*547);
 
 SI3 = u6[i][j-2]*(547*u6[i][j-2]-2522*u6[i][j] + 1922*u6[i][j+2]- 494*u6[i][j+4]) + u6[i][j]*(3443*u6[i][j]-5966*u6[i][j+2]+1602*u6[i][j+4]) + u6[i][j+2]*(2843*u6[i][j+2]-1642*u6[i][j+4]) + u6[i][j+4]*(267*u6[i][j+4]);
 
 SI4 = u6[i][j]*(2107*u6[i][j]-9402*u6[i][j+2] + 7042*u6[i][j+4]- 1854*u6[i][j+6]) + u6[i][j+2]*(11003*u6[i][j+2]-17246*u6[i][j+4]+4642*u6[i][j+6]) + u6[i][j+4]*(7043*u6[i][j+4]-3882*u6[i][j+6]) + u6[i][j+6]*(547*u6[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u6[i][j-6] + (cr5)*u6[i][j-4] - (cr6)*u6[i][j-2]+ (cr7)*u6[i][j]; 
 	
 S2 = (cr1)*u6[i][j-4] - (cr3)*u6[i][j-2] + (cr5)*u6[i][j] + (cr2)*u6[i][j+2]; 
 
 S3 =  (-cr1)*u6[i][j-2] + cr4*u6[i][j] + cr4*u6[i][j+2]- cr1*u6[i][j+4]; 
 
 S4 =  (-cr3)*u6[i][j+4] + cr2*u6[i][j] + cr5*u6[i][j+2]+ cr1*u6[i][j+6]; 
 	
 	
 S_4 = -cr2*u6[i][j+6] + cr5*u6[i][j+4] - cr6*u6[i][j+2]+ cr7*u6[i][j];
 
 S_3 =  cr1*u6[i][j+4] - cr3*u6[i][j+2] + cr5*u6[i][j] + cr2*u6[i][j-2]; 
 
 S_2 = -cr1*u6[i][j+2] + cr4*u6[i][j] + cr4*u6[i][j-2] - cr1*u6[i][j-4] ; 
 
 S_1 =  cr2*u6[i][j]+ cr5*u6[i][j-2] - cr3*u6[i][j-4]+ cr1*u6[i][j-6];
 

 
  bx_L[i][j+1] = (S1*w1)+(S2*w2)+(S3*w3)+ S4*w4;	
 	
  bx_R[i][j-1] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
 	
    
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
   SI1 = u6[i-6][j]*(547*u6[i-6][j]-3882*u6[i-4][j]+4642*u6[i-2][j]-1854*u6[i][j]) + u6[i-4][j]*(7043*u6[i-4][j]-17246*u6[i-2][j]+7042*u6[i][j])+ u6[i-2][j]*(11003*u6[i-2][j]-9402*u6[i][j])+2107*u6[i][j]*u6[i][j];
 
 SI2 = u6[i-4][j]*(267*u6[i-4][j]-1642*u6[i-2][j] + 1602*u6[i][j]-494*u6[i+2][j])+ u6[i-2][j]*(2843*u6[i-2][j]-5966*u6[i][j] + 1922*u6[i+2][j]) +u6[i][j]*(3443*u6[i][j]-2522*u6[i+2][j])+ (u6[i+2][j]*u6[i+2][j]*547);
 
 SI3 = u6[i-2][j]*(547*u6[i-2][j]-2522*u6[i][j] + 1922*u6[i+2][j]- 494*u6[i+4][j]) + u6[i][j]*(3443*u6[i][j]-5966*u6[i+2][j]+1602*u6[i+4][j]) + u6[i+2][j]*(2843*u6[i+2][j]-1642*u6[i+4][j]) + u6[i+4][j]*(267*u6[i+4][j]);
 
 SI4 = u6[i][j]*(2107*u6[i][j]-9402*u6[i+2][j] + 7042*u6[i+4][j]- 1854*u6[i+6][j]) + u6[i+2][j]*(11003*u6[i+2][j]-17246*u6[i+4][j]+4642*u6[i+6][j]) + u6[i+4][j]*(7043*u6[i+4][j]-3882*u6[i+6][j]) + u6[i+6][j]*(547*u6[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u6[i-6][j] + (cr5)*u6[i-4][j] - (cr6)*u6[i-2][j]+ (cr7)*u6[i][j]; 
 	
 S2 = (cr1)*u6[i-4][j] - (cr3)*u6[i-2][j] + (cr5)*u6[i][j] + (cr2)*u6[i+2][j]; 
 
 S3 =  (-cr1)*u6[i-2][j] + cr4*u6[i][j] + cr4*u6[i+2][j]- cr1*u6[i+4][j]; 
 
 S4 =  (-cr3)*u6[i+4][j] + cr2*u6[i][j] + cr5*u6[i+2][j]+ cr1*u6[i+6][j]; 
 	
 	
 S_4 = -cr2*u6[i+6][j] + cr5*u6[i+4][j] - cr6*u6[i+2][j]+ cr7*u6[i][j];
 
 S_3 =  cr1*u6[i+4][j] - cr3*u6[i+2][j] + cr5*u6[i][j] + cr2*u6[i-2][j]; 
 
 S_2 = -cr1*u6[i+2][j] + cr4*u6[i][j] + cr4*u6[i-2][j] - cr1*u6[i-4][j] ; 
 
 S_1 =  cr2*u6[i][j]+ cr5*u6[i-2][j] - cr3*u6[i-4][j]+ cr1*u6[i-6][j];


  BX_L[i+1][j] = (S1*w1)+(S2*w2)+(S3*w3)+S4*w4;	
 	
 BX_R[i-1][j] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
 	
 	 }

}
return;
}

void WENO8(double u7[][points])
{

for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
{
    
 // if(i==10)
  //cin>>j;
 
 SI1 = u7[i][j-6]*(547*u7[i][j-6]-3882*u7[i][j-4]+4642*u7[i][j-2]-1854*u7[i][j]) + u7[i][j-4]*(7043*u7[i][j-4]-17246*u7[i][j-2]+7042*u7[i][j])+ u7[i][j-2]*(11003*u7[i][j-2]-9402*u7[i][j])+2107*u7[i][j]*u7[i][j];
 
 SI2 = u7[i][j-4]*(267*u7[i][j-4]-1642*u7[i][j-2] + 1602*u7[i][j]-494*u7[i][j+2])+ u7[i][j-2]*(2843*u7[i][j-2]-5966*u7[i][j] + 1922*u7[i][j+2]) +u7[i][j]*(3443*u7[i][j]-2522*u7[i][j+2])+ (u7[i][j+2]*u7[i][j+2]*547);
 
 SI3 = u7[i][j-2]*(547*u7[i][j-2]-2522*u7[i][j] + 1922*u7[i][j+2]- 494*u7[i][j+4]) + u7[i][j]*(3443*u7[i][j]-5966*u7[i][j+2]+1602*u7[i][j+4]) + u7[i][j+2]*(2843*u7[i][j+2]-1642*u7[i][j+4]) + u7[i][j+4]*(267*u7[i][j+4]);
 
 SI4 = u7[i][j]*(2107*u7[i][j]-9402*u7[i][j+2] + 7042*u7[i][j+4]- 1854*u7[i][j+6]) + u7[i][j+2]*(11003*u7[i][j+2]-17246*u7[i][j+4]+4642*u7[i][j+6]) + u7[i][j+4]*(7043*u7[i][j+4]-3882*u7[i][j+6]) + u7[i][j+6]*(547*u7[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u7[i][j-6] + (cr5)*u7[i][j-4] - (cr6)*u7[i][j-2]+ (cr7)*u7[i][j]; 
 	
 S2 = (cr1)*u7[i][j-4] - (cr3)*u7[i][j-2] + (cr5)*u7[i][j] + (cr2)*u7[i][j+2]; 
 
 S3 =  (-cr1)*u7[i][j-2] + cr4*u7[i][j] + cr4*u7[i][j+2]- cr1*u7[i][j+4]; 
 
 S4 =  (-cr3)*u7[i][j+4] + cr2*u7[i][j] + cr5*u7[i][j+2]+ cr1*u7[i][j+6]; 
 	
 	
 S_4 = -cr2*u7[i][j+6] + cr5*u7[i][j+4] - cr6*u7[i][j+2]+ cr7*u7[i][j];
 
 S_3 =  cr1*u7[i][j+4] - cr3*u7[i][j+2] + cr5*u7[i][j] + cr2*u7[i][j-2]; 
 
 S_2 = -cr1*u7[i][j+2] + cr4*u7[i][j] + cr4*u7[i][j-2] - cr1*u7[i][j-4] ; 
 
 S_1 =  cr2*u7[i][j]+ cr5*u7[i][j-2] - cr3*u7[i][j-4]+ cr1*u7[i][j-6];
 

 
                                                           
  velocity_zL[i][j+1] = S1*w1+S2*w2+S3*w3+ S4*w4;	
 	
  velocity_zR[i][j-1] = S_1*w_1+S_2*w_2+S_3*w_3+ S_4*w_4;	
 	


 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = u7[i-6][j]*(547*u7[i-6][j]-3882*u7[i-4][j]+4642*u7[i-2][j]-1854*u7[i][j]) + u7[i-4][j]*(7043*u7[i-4][j]-17246*u7[i-2][j]+7042*u7[i][j])+ u7[i-2][j]*(11003*u7[i-2][j]-9402*u7[i][j])+2107*u7[i][j]*u7[i][j];
 
 SI2 = u7[i-4][j]*(267*u7[i-4][j]-1642*u7[i-2][j] + 1602*u7[i][j]-494*u7[i+2][j])+ u7[i-2][j]*(2843*u7[i-2][j]-5966*u7[i][j] + 1922*u7[i+2][j]) +u7[i][j]*(3443*u7[i][j]-2522*u7[i+2][j])+ (u7[i+2][j]*u7[i+2][j]*547);
 
 SI3 = u7[i-2][j]*(547*u7[i-2][j]-2522*u7[i][j] + 1922*u7[i+2][j]- 494*u7[i+4][j]) + u7[i][j]*(3443*u7[i][j]-5966*u7[i+2][j]+1602*u7[i+4][j]) + u7[i+2][j]*(2843*u7[i+2][j]-1642*u7[i+4][j]) + u7[i+4][j]*(267*u7[i+4][j]);
 
 SI4 = u7[i][j]*(2107*u7[i][j]-9402*u7[i+2][j] + 7042*u7[i+4][j]- 1854*u7[i+6][j]) + u7[i+2][j]*(11003*u7[i+2][j]-17246*u7[i+4][j]+4642*u7[i+6][j]) + u7[i+4][j]*(7043*u7[i+4][j]-3882*u7[i+6][j]) + u7[i+6][j]*(547*u7[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u7[i-6][j] + (cr5)*u7[i-4][j] - (cr6)*u7[i-2][j]+ (cr7)*u7[i][j]; 
 	
 S2 = (cr1)*u7[i-4][j] - (cr3)*u7[i-2][j] + (cr5)*u7[i][j] + (cr2)*u7[i+2][j]; 
 
 S3 =  (-cr1)*u7[i-2][j] + cr4*u7[i][j] + cr4*u7[i+2][j]- cr1*u7[i+4][j]; 
 
 S4 =  (-cr3)*u7[i+4][j] + cr2*u7[i][j] + cr5*u7[i+2][j]+ cr1*u7[i+6][j]; 
 	
 	
 S_4 = -cr2*u7[i+6][j] + cr5*u7[i+4][j] - cr6*u7[i+2][j]+ cr7*u7[i][j];
 
 S_3 =  cr1*u7[i+4][j] - cr3*u7[i+2][j] + cr5*u7[i][j] + cr2*u7[i-2][j]; 
 
 S_2 = -cr1*u7[i+2][j] + cr4*u7[i][j] + cr4*u7[i-2][j] - cr1*u7[i-4][j] ; 
 
 S_1 =  cr2*u7[i][j]+ cr5*u7[i-2][j] - cr3*u7[i-4][j]+ cr1*u7[i-6][j];


v_zL[i+1][j] = S1*w1+S2*w2+S3*w3+ S4*w4;	
 	
  v_zR[i-1][j] = S_1*w_1+S_2*w_2+S_3*w_3+ S_4*w_4;	
 
 
 	 }
 	// return;
}
}

void WENO9(double u8[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=6; i<=size-6; i=i+2)
 {
 	
for(j=6; j<=size-6; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 


SI1 = u8[i][j-6]*(547*u8[i][j-6]-3882*u8[i][j-4]+4642*u8[i][j-2]-1854*u8[i][j]) + u8[i][j-4]*(7043*u8[i][j-4]-17246*u8[i][j-2]+7042*u8[i][j])+ u8[i][j-2]*(11003*u8[i][j-2]-9402*u8[i][j])+2107*u8[i][j]*u8[i][j];
 
 SI2 = u8[i][j-4]*(267*u8[i][j-4]-1642*u8[i][j-2] + 1602*u8[i][j]-494*u8[i][j+2])+ u8[i][j-2]*(2843*u8[i][j-2]-5966*u8[i][j] + 1922*u8[i][j+2]) +u8[i][j]*(3443*u8[i][j]-2522*u8[i][j+2])+ (u8[i][j+2]*u8[i][j+2]*547);
 
 SI3 = u8[i][j-2]*(547*u8[i][j-2]-2522*u8[i][j] + 1922*u8[i][j+2]- 494*u8[i][j+4]) + u8[i][j]*(3443*u8[i][j]-5966*u8[i][j+2]+1602*u8[i][j+4]) + u8[i][j+2]*(2843*u8[i][j+2]-1642*u8[i][j+4]) + u8[i][j+4]*(267*u8[i][j+4]);
 
 SI4 = u8[i][j]*(2107*u8[i][j]-9402*u8[i][j+2] + 7042*u8[i][j+4]- 1854*u8[i][j+6]) + u8[i][j+2]*(11003*u8[i][j+2]-17246*u8[i][j+4]+4642*u8[i][j+6]) + u8[i][j+4]*(7043*u8[i][j+4]-3882*u8[i][j+6]) + u8[i][j+6]*(547*u8[i][j+6]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u8[i][j-6] + (cr5)*u8[i][j-4] - (cr6)*u8[i][j-2]+ (cr7)*u8[i][j]; 
 	
 S2 = (cr1)*u8[i][j-4] - (cr3)*u8[i][j-2] + (cr5)*u8[i][j] + (cr2)*u8[i][j+2]; 
 
 S3 =  (-cr1)*u8[i][j-2] + (cr4)*u8[i][j] + (cr4)*u8[i][j+2]- (cr1)*u8[i][j+4]; 
 
 S4 =  (-cr3)*u8[i][j+4] + (cr2)*u8[i][j] + (cr5)*u8[i][j+2]+ (cr1)*u8[i][j+6]; 
 	
 	
 S_4 = (-cr2)*u8[i][j+6] + (cr5)*u8[i][j+4] - (cr6)*u8[i][j+2]+ (cr7)*u8[i][j];
 
 S_3 =  cr1*u8[i][j+4] - cr3*u8[i][j+2] + cr5*u8[i][j] + (cr2)*u8[i][j-2]; 
 
 S_2 = -cr1*u8[i][j+2] + cr4*u8[i][j] + cr4*u8[i][j-2] - (cr1)*u8[i][j-4] ; 
 
 S_1 =  cr2*u8[i][j]+ cr5*u8[i][j-2] - cr3*u8[i][j-4]+ cr1*u8[i][j-6];
 
 
  bz_L[i][j+1] =  (S1*w1)+(S2*w2)+(S3*w3)+ S4*w4;	
 	
  bz_R[i][j-1] =  (S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
 	
    
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
   SI1 = u8[i-6][j]*(547*u8[i-6][j]-3882*u8[i-4][j]+4642*u8[i-2][j]-1854*u8[i][j]) + u8[i-4][j]*(7043*u8[i-4][j]-17246*u8[i-2][j]+7042*u8[i][j])+ u8[i-2][j]*(11003*u8[i-2][j]-9402*u8[i][j])+2107*u8[i][j]*u8[i][j];
 
 SI2 = u8[i-4][j]*(267*u8[i-4][j]-1642*u8[i-2][j] + 1602*u8[i][j]-494*u8[i+2][j])+ u8[i-2][j]*(2843*u8[i-2][j]-5966*u8[i][j] + 1922*u8[i+2][j]) +u8[i][j]*(3443*u8[i][j]-2522*u8[i+2][j])+ (u8[i+2][j]*u8[i+2][j]*547);
 
 SI3 = u8[i-2][j]*(547*u8[i-2][j]-2522*u8[i][j] + 1922*u8[i+2][j]- 494*u8[i+4][j]) + u8[i][j]*(3443*u8[i][j]-5966*u8[i+2][j]+1602*u8[i+4][j]) + u8[i+2][j]*(2843*u8[i+2][j]-1642*u8[i+4][j]) + u8[i+4][j]*(267*u8[i+4][j]);
 
 SI4 = u8[i][j]*(2107*u8[i][j]-9402*u8[i+2][j] + 7042*u8[i+4][j]- 1854*u8[i+6][j]) + u8[i+2][j]*(11003*u8[i+2][j]-17246*u8[i+4][j]+4642*u8[i+6][j]) + u8[i+4][j]*(7043*u8[i+4][j]-3882*u8[i+6][j]) + u8[i+6][j]*(547*u8[i+6][j]);
 
 
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 omega4= gamma4/((elipson + SI4)*(elipson + SI4));
 
 
 omega_1= gamma4/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma3/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma2/((elipson+SI3)*(elipson+SI3));
 omega_4= gamma1/((elipson+SI4)*(elipson+SI4));
 
 
 w1=omega1/(omega1+omega2+omega3+omega4);
 w2=omega2/(omega1+omega2+omega3+omega4);
 w3=omega3/(omega1+omega2+omega3+omega4);
 w4=omega4/(omega1+omega2+omega3+omega4);
 
 w_1=omega_1/(omega_1+omega_2+omega_3+omega_4);
 w_2=omega_2/(omega_1+omega_2+omega_3+omega_4);
 w_3=omega_3/(omega_1+omega_2+omega_3+omega_4);
 w_4=omega_4/(omega_1+omega_2+omega_3+omega_4);
 
 
 
 S1 = (-cr2)*u8[i-6][j] + (cr5)*u8[i-4][j] - (cr6)*u8[i-2][j]+ (cr7)*u8[i][j]; 
 	
 S2 = (cr1)*u8[i-4][j] - (cr3)*u8[i-2][j] + (cr5)*u8[i][j] + (cr2)*u8[i+2][j]; 
 
 S3 =  (-cr1)*u8[i-2][j] + cr4*u8[i][j] + cr4*u8[i+2][j]- cr1*u8[i+4][j]; 
 
 S4 =  (-cr3)*u8[i+4][j] + cr2*u8[i][j] + cr5*u8[i+2][j]+ cr1*u8[i+6][j]; 
 	
 	
 S_4 = -cr2*u8[i+6][j] + cr5*u8[i+4][j] - cr6*u8[i+2][j]+ cr7*u8[i][j];
 
 S_3 =  cr1*u8[i+4][j] - cr3*u8[i+2][j] + cr5*u8[i][j] + cr2*u8[i-2][j]; 
 
 S_2 = -cr1*u8[i+2][j] + cr4*u8[i][j] + cr4*u8[i-2][j] - cr1*u8[i-4][j] ; 
 
 S_1 =  cr2*u8[i][j]+ cr5*u8[i-2][j] - cr3*u8[i-4][j]+ cr1*u8[i-6][j];

 
  BZ_L[i+1][j] =  (S1*w1)+(S2*w2)+(S3*w3)+ S4*w4;	
 	
 BZ_R[i-1][j] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3)+ S_4*w_4;	
 
 	
 /*	  bz_L[i][j+1]  =    u8[i][j] + (u8[i][j+2]-u8[i][j])*0.5*VANLEER1(u8[i][j],u8[i][j+2],u8[i][j-2]);	
 	
   bz_R[i][j-1]  =    u8[i][j] - (u8[i][j+2]-u8[i][j])*0.5*VANLEER1(u8[i][j],u8[i][j+2],u8[i][j-2]);

	
   BZ_L[i][j]  =    u8[i][j] + (u8[i+2][j]-u8[i][j])*0.5*VANLEER1(u8[i][j],u8[i+2][j],u8[i-2][j]);	
 	
   BZ_R[i][j]  =    u8[i][j] - (u8[i+2][j]-u8[i][j])*0.5*VANLEER1(u8[i][j],u8[i+2][j],u8[i-2][j]);
*/	 //cout<<BZ_L[i+1][j]<<endl;
	 
	  }

}
return;
}

