//----Mag code-----
//--x=i,y=j,z=k---------
////---X[i][j][k][field][time]---
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <omp.h>
#include <complex>
using namespace std;

const double delta_x =0.1;
const int n_max= 100;
double delta_t =0.002;//0.002
const int t_max=16000;
const int No_field=3;
double E=0.0,KK=0.0,V=0.0;
const int NOTHREADS=6;
double X[n_max][n_max][3][3];
double Y[n_max][n_max][3][3];
double points[6][2];
const double PI = 2.0*acos(0);
double energy=0.0,filenumber=0.0;
double dtdt=0.0,dtdty=0.0,dtpdtpy=0.0,d2pd2p=0.0,d1pd1p=0.0,d2yd2y=0.0,d1yd1y=0.0,dtpdtp=0.0,lambday =0.0,lambda =0.0,d22pd2p=0.0;

double d1p(int i,int j, int f, int t);
double d2p(int i,int j, int f, int t);
double d11p(int i,int j, int f, int t);
double d22p(int i,int j, int f, int t);

double d1y(int i,int j, int f, int t);
double d2y(int i,int j, int f, int t);
double d11y(int i,int j, int f, int t);
double d22y(int i,int j, int f, int t);
double X1=0.0,X2=0.0,X3=0.0,Y1=0.0,Y2=0.0,Y3=0.0;
double X11=0.0,X22=0.0,X33=0.0,Y11=0.0,Y22=0.0,Y33=0.0;
double F(double r);

void ConfigEnergy();
void ConfigEnergyy();
const double ci =(n_max/2.0)*delta_x;
const double cj =(n_max/2.0)*delta_x;
double x=0.0,y=0.0,r=0.0,theta=0.0;
double xa=0.0,ya=0.0,ra=0.0,thetaa=0.0;
const double gam=0.5*PI,M=1,gamy=-0.5*PI,sepx=0.25;
const double Coupling=2.0;
double e=0.0;
double alpha[3];
double alphay[3];
void filename(int i, char name2[], char file[]);
double dum =0.0;
double ap=0.0,apy=0.0;
char file[]="xx_0000";

const double mu=-0.9;
double a1=0.0;

const int new_time = 2000000000;
const double m=0.5,D=1.0,Dy=-1.0;
double DTS=0.0,DTS1=0.0,DTS2=0.0,DTS3=0.0;
double DTSY=0.0,DTS1Y=0.0,DTS2Y=0.0,DTS3Y=0.0;
int e1=0,e2=0,e3=0;
int main( )
{
  //-------read me file--------------
ofstream ReadMe;
ReadMe.open("Read-Me");
 ReadMe << "2D Magnetic Skyrmion gradient flow simulation  \n" << " delta x = " <<delta_x<<"\n" <<"delta_t = " << delta_t<<"\n"<<" n_max = " <<n_max<<"\n"<<"\n delta_t = "<<delta_t<<"\n t_max = " <<t_max<< "\n "; 

 ReadMe.close();
 //------read me end -----------------
 
 ofstream Hopfenergy;
Hopfenergy.open("Hopfenergy");
 complex<double> I(0,1.0);
 //--Set fields to zero
for(int i =0; i < n_max; ++i)
  {
 for( int j = 0; j < n_max; ++j)
   {
   X[i][j][0][0]=X[i][j][0][1]=X[i][j][0][2] = 0.0;
   X[i][j][1][0]=X[i][j][1][1]=X[i][j][1][2] = 0.0;
   X[i][j][2][0]=X[i][j][2][1]=X[i][j][2][2] = 1.0;
 
   Y[i][j][0][1] = 0.0;
   Y[i][j][1][1] = 0.0;
   Y[i][j][2][1] = 1.0;
	
 }
 }
//--positions of Skyrmions
 points[0][0]=-3.5;
 points[0][1]=0.0;
points[1][0]=-2.5;
 points[1][1]=3.0;
points[2][0]=-2.5;
 points[2][1]=-3.0;

points[3][0]=+1.5;
 points[3][1]=0.0;
points[4][0]=+2.5;
 points[4][1]=3.0;
points[5][0]=+2.5;
 points[5][1]=-3.0;
 
//----------------- inital configuration  ----------------
 for(int k=0;k<6;++k)
   {
for(int i =1; i < n_max-1; ++i)
  {x=i*delta_x-ci+sepx+points[k][0],xa=i*delta_x-ci-sepx+points[k][0];
 for( int j = 1; j < n_max-1; ++j)
  {y=j*delta_x-cj+points[k][1]+0.1,ya=j*delta_x-cj+points[k][1]-0.1;
   r=sqrt(pow(x,2)+pow(y,2));
   theta=atan2(y,x);
   ra=sqrt(pow(xa,2)+pow(ya,2));
   thetaa=atan2(ya,xa);

   X1 = sin(F(r))*cos(M*theta+gam);
   X2 = sin(F(r))*sin(M*theta+gam);
   X3 = cos(F(r));
 
   Y1 = sin(F(ra))*cos(M*thetaa+gamy);
   Y2 = sin(F(ra))*sin(M*thetaa+gamy);
   Y3 = cos(F(ra));
	
   complex<double> IXC(X[i][j][0][1]/(1.0+X[i][j][2][1]),X[i][j][1][1]/(1.0+X[i][j][2][1])); //-initial values
   complex<double> IYC(Y[i][j][0][1]/(1.0+Y[i][j][2][1]),Y[i][j][1][1]/(1.0+Y[i][j][2][1]));
   
   complex<double> XC(X1/(1.0+X3),X2/(1.0+X3));
   complex<double> YC(Y1/(1.0+Y3),Y2/(1.0+Y3));

   complex<double> WX =XC+IXC;
   complex<double> WY =YC+IYC;
   
   X11=real((WX+conj(WX))/(1.0+WX*conj(WX)));
   X22=real(I*(conj(WX)-WX)/(1.0+WX*conj(WX)));
   X33=real((1.0-conj(WX)*WX)/(1.0+WX*conj(WX)));

   Y11=real((WY+conj(WY))/(1.0+WY*conj(WY)));
   Y22=real(I*(conj(WY)-WY)/(1.0+WY*conj(WY)));
   Y33=real((1.0-conj(WY)*WY)/(1.0+WY*conj(WY)));

  
   X[i][j][0][1]=X11;
   X[i][j][1][1]=X22;
   X[i][j][2][1]=X33;

   Y[i][j][0][1]=Y11;
   Y[i][j][1][1]=Y22;
   Y[i][j][2][1]=Y33;
   
   
 }
 }
   }

ConfigEnergy();
	 
//--Update all fields--
for(int i =1; i < n_max-1; ++i)
  {
 for( int j = 1; j < n_max-1; ++j)
   {
   
   X[i][j][0][2]= X[i][j][0][0]=X[i][j][0][1]; 
   X[i][j][1][2]= X[i][j][1][0]=X[i][j][1][1]; 
   X[i][j][2][2]= X[i][j][2][0]=X[i][j][2][1]; 
 
   Y[i][j][0][2]= Y[i][j][0][0]=Y[i][j][0][1];
   Y[i][j][1][2]= Y[i][j][1][0]=Y[i][j][1][1];
   Y[i][j][2][2]= Y[i][j][2][0]=Y[i][j][2][1];
	
 
 }
 }
ConfigEnergy();
	 cout << "X2 Field \n ";
cout << "Energy = " <<E<<", kintetic =" <<KK<<", Both ="<<E+KK<<", delta_t  =" << delta_t<<endl;
 //-----  end of initial data ------------------


 //-----------------Data from previouse 2d ----------------

 ifstream dat("./f1_0006");
 for(int i =1+0*15; i < n_max-1-0*15; ++i)
 {
 for( int j = 1+0*15; j < n_max-1-0*15; ++j)
 {
 dat >>a1>> X[i][j][0][0]>>X[i][j][1][0] >>X[i][j][2][0]>> Y[i][j][0][0]>>Y[i][j][1][0] >>Y[i][j][2][0];
 
X[i][j][0][1]=X[i][j][0][0];
 X[i][j][1][1]=X[i][j][1][0];
 X[i][j][2][1]=X[i][j][2][0]; 
 
 X[i][j][0][2]=X[i][j][0][0];
 X[i][j][1][2]=X[i][j][1][0];
 X[i][j][2][2]=X[i][j][2][0]; 

 Y[i][j][0][1]=Y[i][j][0][0];
 Y[i][j][1][1]=Y[i][j][1][0];
 Y[i][j][2][1]=Y[i][j][2][0]; 
 
 Y[i][j][0][2]=Y[i][j][0][0];
 Y[i][j][1][2]=Y[i][j][1][0];
 Y[i][j][2][2]=Y[i][j][2][0];
 
 }
 }	 
dat.close();


  
 //------------time evolution -----------------------------------
 for(int t=0; t<t_max;++t)
   { 
      
     for(int i=1+0*25; i < n_max-1-0*25;++i)
       {
#pragma omp parallel num_threads(NOTHREADS)
	       {
		  #pragma omp for nowait	\
		    private(lambda,alpha,d2pd2p,d1pd1p,dtpdtp,ap,d2yd2y=0.0,d1yd1y=0.0,dtpdtpy) \
  schedule(guided)
      for( int j = 1+0*35;j < n_max-1-0*35; ++j)
		 {
	
// --- calculating skyrme dot products ------------------
		   // ---reset all dot's ------
		 d2pd2p=0.0,d1pd1p=0.0,dtpdtp=0.0,lambda =0.0;
		 d2yd2y=0.0,d1yd1y=0.0,dtpdtpy=0.0,lambday =0.0;
		
		   //-----------dot products ----------
  for(int f =0; f <No_field; ++f)
		     {
	  d2pd2p  += pow(d2p(i,j,f,1),2);
	  
	  d1pd1p  += pow(d1p(i,j,f,1),2);

	  d2yd2y  += pow(d2y(i,j,f,1),2);
	  
	  d1yd1y  += pow(d1y(i,j,f,1),2);
	  
	  dtpdtp +=pow(((X[i][j][f][1]-X[i][j][f][0])/delta_t),2);

	  dtpdtpy+=pow(((Y[i][j][f][1]-Y[i][j][f][0])/delta_t),2);
		     }
 	   //------------calc next time value ---------------------- 
  ap=0.0;
  apy=0.0;  
  for(int q=0; q<No_field; ++q)
    {if(q==No_field-1){e=1.0;}else{e=0;}
if(q==0){e1=1.0;}else{e1=0;}
if(q==1){e2=1.0;}else{e2=0;}
if(q==2){e3=1.0;}else{e3=0;}

 alpha[q] = pow(m,2)*e+Coupling*Y[i][j][q][1]
   -2.0*D*((d2p(i,j,2,1))*e1+(-d1p(i,j,2,1))*e2+(d1p(i,j,1,1)-d2p(i,j,0,1))*e3);
  
  ap+=X[i][j][q][1]*alpha[q];

alphay[q] = pow(m,2)*e+Coupling*X[i][j][q][1]
   -2.0*Dy*((d2y(i,j,2,1))*e1+(-d1y(i,j,2,1))*e2+(d1y(i,j,1,1)-d2y(i,j,0,1))*e3);
  
  apy+=Y[i][j][q][1]*alphay[q];
		     } 
 
		       	       
		       //----Lagrange multiplyer ---------------
 lambda = -dtpdtp +(d1pd1p+d2pd2p)-ap; 
 lambday=-dtpdtpy +(d1yd1y+d2yd2y)-apy;
		   //---------next time step --------------
	for(int s=0; s<No_field; ++s)
		     {
X[i][j][s][2]= ((d11p(i,j,s,1)+d22p(i,j,s,1))+alpha[s] +X[i][j][s][1]*lambda )*pow(delta_t,2) +2.0*X[i][j][s][1]-X[i][j][s][0]; 
Y[i][j][s][2]= ((d11y(i,j,s,1)+d22y(i,j,s,1))+alphay[s]+Y[i][j][s][1]*lambday)*pow(delta_t,2) +2.0*Y[i][j][s][1]-Y[i][j][s][0]; 
		
             }


		       	 }}
 
       }
       
       //-----Energy -------------
      if(t%20 ==0)
	 {
	 ConfigEnergy();
	 cout << "X Field \n ";
cout << "Energy = " <<E<<", kintetic =" <<KK<<", Both ="<<E+KK<<", delta_t  =" << delta_t<<endl;  
 Hopfenergy <<t<<"    "<<E<<"    "<<KK<<"    " <<E+KK<<"\n";

ConfigEnergyy();
 cout << "Y Field \n ";
cout << "Energy = " <<E<<", kintetic =" <<KK<<", Both ="<<E+KK<<", delta_t  =" << delta_t<<endl;  
 Hopfenergy <<t<<"    "<<E<<"    "<<KK<<"    " <<E+KK<<"\n";

	 }
//------regular field output.
       if( (t/(t_max/5.0))-int(t/int(t_max/5)) ==0)
	 {
	   E=0.0,KK=0.0;
	   
	    filenumber++;
	   	   filename(int(filenumber), "f1",file);
	   	    ofstream outf1(file);

for(int i =1; i < n_max-1; ++i)
       {
	 for(int j =1; j<n_max-1; ++j)
	   {
		 outf1 << energy  <<"   "<< X[i][j][0][0]<<"     "<<X[i][j][1][0] <<"    " << X[i][j][2][0]<<"   "<< Y[i][j][0][0]<<"     "<<Y[i][j][1][0] <<"    " << Y[i][j][2][0] <<"\n";     
	   }
	       }}
    //outf1.close();
      
 //----update field ------------------------------------
 #pragma omp parallel num_threads(NOTHREADS)
	       {
		  #pragma omp for nowait	\
		    private(DTS,DTS1)		\
      schedule(guided)
 for(int i=1+0*25; i < n_max-1-0*25; ++i)
       {
	 for(int j =1+0*25;j < n_max-1-0*25; ++j)
	   {
	     
   DTS = sqrt(pow(X[i][j][0][2],2)+pow(X[i][j][1][2],2)+pow(X[i][j][2][2],2));
   DTSY =sqrt(pow(Y[i][j][0][2],2)+pow(Y[i][j][1][2],2)+pow(Y[i][j][2][2],2));
   if((t%20 ==0 && t > 50) || (t%10 ==0 && t < 50))
         {
				
	   	X[i][j][0][1]=X[i][j][0][0]=(X[i][j][0][2])/DTS;
	   	X[i][j][1][1]=X[i][j][1][0]=(X[i][j][1][2])/DTS;
		X[i][j][2][1]=X[i][j][2][0]=(X[i][j][2][2])/DTS;

		Y[i][j][0][1]=Y[i][j][0][0]=(Y[i][j][0][2])/DTSY;
	   	Y[i][j][1][1]=Y[i][j][1][0]=(Y[i][j][1][2])/DTSY;
		Y[i][j][2][1]=Y[i][j][2][0]=(Y[i][j][2][2])/DTSY;
			 
				} 
			     else
			       
			       {

		   	     DTS1 = sqrt(pow(X[i][j][0][1],2)+pow(X[i][j][1][1],2)+pow(X[i][j][2][1],2));
			     DTS1Y= sqrt(pow(Y[i][j][0][1],2)+pow(Y[i][j][1][1],2)+pow(Y[i][j][2][1],2));

			     X[i][j][0][0]=(X[i][j][0][1])/DTS1;
			     X[i][j][0][1]=(X[i][j][0][2])/DTS;
			     X[i][j][1][0]=(X[i][j][1][1])/DTS1;
			     X[i][j][1][1]=(X[i][j][1][2])/DTS;
			     X[i][j][2][0]=(X[i][j][2][1])/DTS1;
			     X[i][j][2][1]=(X[i][j][2][2])/DTS;
					     
			     Y[i][j][0][0]=(Y[i][j][0][1])/DTS1Y;
			     Y[i][j][0][1]=(Y[i][j][0][2])/DTSY;
			     Y[i][j][1][0]=(Y[i][j][1][1])/DTS1Y;
			     Y[i][j][1][1]=(Y[i][j][1][2])/DTSY;
			     Y[i][j][2][0]=(Y[i][j][2][1])/DTS1Y;
			     Y[i][j][2][1]=(Y[i][j][2][2])/DTSY;
			     
			     }

	       	       }}
  }

   }
 

 
 //-------final energy/configuration -----------------------------
   E=0.0,KK=0.0;
  
    filenumber++;
 	   filename(int(filenumber), "f1",file);
 	    ofstream outf1(file);
	   
for(int i =1; i < n_max-1; ++i)
       {x=i*delta_x-ci;
	 for(int j =1; j<n_max-1; ++j)
	   {y=j*delta_x-ci;	    

		 outf1 <<  acos(X[i][j][0][0]*sqrt(1.0-pow(mu,2))+mu*X[i][j][2][0] )/acos(-1.0*mu)  <<"   "<< X[i][j][0][0]<<"     "<<X[i][j][1][0] <<"    " << X[i][j][2][0]<<"   "<< Y[i][j][0][0]<<"     "<<Y[i][j][1][0] <<"    " << Y[i][j][2][0] <<"\n";
	    
	   }
	   
       }
ofstream SkyrmionConfig;
SkyrmionConfig.open("SkyrmionConfig");
 double scale=0.4;
for(int i =0; i < n_max; i+=4)
       {x=i*delta_x-ci;
	 for(int j =0; j<n_max; j+=4)
	   {y=j*delta_x-cj;
	    r=sqrt(pow(x,2)+pow(y,2));
	    SkyrmionConfig<<x<<"  "<<y<<"  "<<"0  "<< scale*X[i][j][0][1] <<"  "<<scale*X[i][j][1][1]<<"   "<<scale*X[i][j][2][1]<<"  "<< scale*Y[i][j][0][1] <<"  "<<scale*Y[i][j][1][1]<<"   "<<scale*Y[i][j][2][1]<<"\n";
	   }
	   
       }
   //outf1.close();
 Hopfenergy.close();
SkyrmionConfig.close();
 cout<< " Final energy = "<<E<"\n";
 


 
 return 0;
	   }

void ConfigEnergy()
{ E=0.0,KK=0.0;
for(int i =1; i < n_max-1; ++i)
       {
	 for(int j =1; j<n_max-1; ++j)
	   {
	     d2pd2p=0.0,d1pd1p=0.0,dtdt=0.0;
	     for(int f =0; f <3; ++f)
	       {
		 dtdt +=pow(((X[i][j][f][2]-X[i][j][f][0])/(2.0*delta_t)),2);
		 d1pd1p  +=pow(d1p(i,j,f,0),2);
				
		 d2pd2p  +=pow(d2p(i,j,f,0),2);	
	       }
	     
	     E+=energy=((0.5*(d1pd1p+d2pd2p)+pow(m,2)*pow(1.0-X[i][j][2][0],1)
 +D*((d2p(i,j,2,0))*X[i][j][0][0]+(-d1p(i,j,2,0))*X[i][j][1][0]+(d1p(i,j,1,0)-d2p(i,j,0,0))*X[i][j][2][0])
 )*pow(delta_x,2)); 
 KK+= dtdt*pow(delta_x,2);
 
	       }
	       }
}

void ConfigEnergyy()
{ E=0.0,KK=0.0;
for(int i =1; i < n_max-1; ++i)
       {
	 for(int j =1; j<n_max-1; ++j)
	   {
	     d2yd2y=0.0,d1yd1y=0.0,dtdty=0.0;
	     for(int f =0; f <3; ++f)
	       {
		 dtdty +=pow(((Y[i][j][f][2]-Y[i][j][f][0])/(2.0*delta_t)),2);
		 d1yd1y  +=pow(d1y(i,j,f,1),2);
				
		 d2yd2y  +=pow(d2y(i,j,f,1),2);
			
	       }
	     E+=energy=((0.5*(d1yd1y+d2yd2y)+pow(m,2)*pow(1.0-Y[i][j][2][0],1)
 +Dy*((d2y(i,j,2,1))*Y[i][j][0][0]+(-d1y(i,j,2,1))*Y[i][j][1][0]+(d1y(i,j,1,1)-d2y(i,j,0,1))*Y[i][j][2][0])
 )*pow(delta_x,2)); 
 KK+= dtdty*pow(delta_x,2);
	       }
	       }
}
 
void filename(int i, char name2[], char file[])
{
  int i1,i2,i3,i4;

  i1=int(floor(i*1.0/1000.0));
  i=i-i1*1000;
  i2=int(floor(i*1.0/100.0));
  i=i-i2*100;
  i3=int(floor(i*1.0/10.0));
  i=i-i3*10;
  i4=i;

  file[3]=char(48+i1);
  file[4]=char(48+i2);
  file[5]=char(48+i3);
  file[6]=char(48+i4);

  file[0]=name2[0];
  file[1]=name2[1];
}

 double F(double r)
 {return PI*exp(-r);
		 // return PI*((0.5*(n_max*delta_x-delta_x)-r)/(0.5*n_max*delta_x-delta_x)); //--Normal mass case
  // double a=1.6,f=0.0;
 
  // {f= PI*((a-r)/a);}
   //  return f;
}

//----- Higher order derivatives ----------------------------
double d1p(int i,int j,int f, int t)
{
  if(i ==1||i==n_max-2||j ==1||j==n_max-2)
    {
      return (X[i+1][j][f][t]-X[i-1][j][f][t])/(2.0*delta_x);
    } 
  else
    {
      return  (8.0*X[i+1][j][f][t]-8.0*X[i-1][j][f][t]-X[i+2][j][f][t]+X[i-2][j][f][t])/(12.0*delta_x);
    }
}

double d2p(int i,int j,int f, int t)
{
  if(i ==1||i==n_max-2||j ==1||j==n_max-2)
    {return (X[i][j+1][f][t]-X[i][j-1][f][t])/(2.0*delta_x);}
  else
    {
      return (8.0*X[i][j+1][f][t]-8.0*X[i][j-1][f][t]-X[i][j+2][f][t]+X[i][j-2][f][t])/(12.0*delta_x);
    }
}


double d11p(int i,int j,int f, int t)
{
 if(i ==1||i==n_max-2||j ==1||j==n_max-2)
    {return (X[i+1][j][f][t]-2.0*X[i][j][f][t]+X[i-1][j][f][t])/(pow(delta_x,2));}
  else
    {
      return (16.0*X[i+1][j][f][t]+16.0*X[i-1][j][f][t]-30.0*X[i][j][f][t]-X[i+2][j][f][t]-X[i-2][j][f][t])/(12.0*pow(delta_x,2));
    }
}

double  d22p(int i,int j,int f, int t)
{
  if(i ==1||i==n_max-2||j ==1||j==n_max-2)
    {return (X[i][j+1][f][t]-2.0*X[i][j][f][t]+X[i][j-1][f][t])/(pow(delta_x,2));}
   else
     {
       return (16.0*X[i][j+1][f][t]+16.0*X[i][j-1][f][t]-30.0*X[i][j][f][t]-X[i][j+2][f][t]-X[i][j-2][f][t])/(12.0*pow(delta_x,2));
     }
}
//----Y field derivatives ---
double d1y(int i,int j,int f, int t)
{
  if(i ==1||i==n_max-2||j ==1||j==n_max-2)
    {
      return (Y[i+1][j][f][t]-Y[i-1][j][f][t])/(2.0*delta_x);
    } 
  else
    {
      return  (8.0*Y[i+1][j][f][t]-8.0*Y[i-1][j][f][t]-Y[i+2][j][f][t]+Y[i-2][j][f][t])/(12.0*delta_x);
    }
}

double d2y(int i,int j,int f, int t)
{
  if(i ==1||i==n_max-2||j ==1||j==n_max-2)
    {return (Y[i][j+1][f][t]-Y[i][j-1][f][t])/(2.0*delta_x);}
  else
    {
      return (8.0*Y[i][j+1][f][t]-8.0*Y[i][j-1][f][t]-Y[i][j+2][f][t]+Y[i][j-2][f][t])/(12.0*delta_x);
    }
}


double d11y(int i,int j,int f, int t)
{
 if(i ==1||i==n_max-2||j ==1||j==n_max-2)
    {return (Y[i+1][j][f][t]-2.0*Y[i][j][f][t]+Y[i-1][j][f][t])/(pow(delta_x,2));}
  else
    {
      return (16.0*Y[i+1][j][f][t]+16.0*Y[i-1][j][f][t]-30.0*Y[i][j][f][t]-Y[i+2][j][f][t]-Y[i-2][j][f][t])/(12.0*pow(delta_x,2));
    }
}

double  d22y(int i,int j,int f, int t)
{
  if(i ==1||i==n_max-2||j ==1||j==n_max-2)
    {return (Y[i][j+1][f][t]-2.0*Y[i][j][f][t]+Y[i][j-1][f][t])/(pow(delta_x,2));}
   else
     {
       return (16.0*Y[i][j+1][f][t]+16.0*Y[i][j-1][f][t]-30.0*Y[i][j][f][t]-Y[i][j+2][f][t]-Y[i][j-2][f][t])/(12.0*pow(delta_x,2));
     }
}

