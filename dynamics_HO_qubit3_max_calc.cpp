#include<iostream>
#include</include/armadillo>
#include <cmath>
#include <complex>
// #include <nlopt.hpp>
#include<fstream>
#include <vector>
#include <cstdlib>
#define QICLIB_DONT_USE_NLOPT
//# include <QIClib>
#include </home/pratha/QIClib/QIClib-1.0.3/include/QIClib>
// # include </home/aparajita/QIClib-1.0.3/include/QIClib>

using namespace std;
using namespace arma;
using namespace qic;

complex<double> im ={0,1};


mat aa(int N)
{
mat a(N, N, fill::zeros);

for(int n=0; n<=N-1; n++)
 {
  for(int m=0; m<=N-1; m++)
    {

      if(m==n+1){
             a(n,m)=sqrt(m);
             }     
     }
   }

return a;
}

mat aadagger(int N)
{

mat adagger(N, N, fill::zeros);

for(int n=0; n<=N-1; n++)
 {
  for(int m=0; m<=N-1; m++)
    {
       if(n==m+1){
              adagger(n,m)=sqrt(n);
             }  
  }
   }

return adagger; 
}

cx_mat rh0(int N)
{
cx_vec eigval; cx_mat eigvec;

mat a0=aa(N); mat adag0=aadagger(N);
mat adaga=adag0*a0;

eig_gen(eigval, eigvec, adaga);
cx_vec psi0=eigvec.col(0);
cx_mat rhoa0=psi0*trans(psi0); 


mat rhob0={{0,0},{0,1}};

cx_mat rho0=tensor(rhoa0,rhob0);

return rho0;
}
//--------------------- functions for calculating gamma of anhamonic oscillator--------
double f(double omega_m, double beta){
    return 1.0 / (exp(beta * omega_m) - 1.0);
    
}

double J(double omega_m, double alpha, double capOmega) {
    return alpha * omega_m * exp(-omega_m / capOmega);
}


double calculateGamma(double omega_m, double alpha, double capOmega, double beta) {
    double gamma;
    if (omega_m > pow(10,-8)) {
        gamma = J(omega_m,alpha,capOmega) * (1.0 + f(omega_m, beta));
    } else {
        gamma = J(abs(omega_m),alpha,capOmega) * f(abs(omega_m), beta);
    }
    return gamma;
}

/*complex<double> energy_passive(cx_mat rhob, cx_vec e_n)
{
 cx_vec r_n; cx_mat ket_r_n; // eigenvalues and eigenvectors of rho
 eig_gen(r_n,ket_r_n, rhob);
 complex<double> E_passive = r_n(1)*e_n(0)+r_n(0)*e_n(1);
 return E_passive;
}*/

cx_mat L(double g1, double g2, double F, int N, cx_mat rho0)
{
double w0=1.06;
double beta = 1.0;
double capOmega = 1000;
double alpha = 1.0;
mat I2 = eye(2,2);

mat paulix={{0,1},{1,0}}; cx_mat pauliy={{0,-im},{im,0}}; mat pauliz={{1,0},{0,-1}};
cx_mat splus=(paulix+im*pauliy)/2.0; cx_mat sminus=(paulix-im*pauliy)/2.0;

mat a=aa(N); mat adag=aadagger(N);

//------------------- commutator-----------------------------------------------------------

cx_mat f1=g1*(tensor(a,splus)+tensor(adag,sminus))+g2*(tensor(a*a,splus)+tensor(adag*adag,sminus))+F*(tensor(a,I2)+tensor(adag,I2));
cx_mat com=-im*(f1*rho0-rho0*f1);

//-----------------dissipator--------------------------------------------------------
// double nb=1/(exp(w0*beta)-1.0);

// double nb=0.0 , gamma=w0;
double gamma1=calculateGamma(w0,alpha,capOmega,beta);
double gamma2=calculateGamma(-w0,alpha,capOmega,beta);

cx_mat di1=gamma1*(tensor(a,I2)*rho0*tensor(adag,I2)-0.5*(tensor(adag,I2)*tensor(a,I2)*rho0+rho0*tensor(adag,I2)*tensor(a,I2)));
cx_mat di2=gamma2*(tensor(adag,I2)*rho0*tensor(a,I2)-0.5*(tensor(a,I2)*tensor(adag,I2)*rho0+rho0*tensor(a,I2)*tensor(adag,I2)));
// std::cout << gamma1 <<"," <<gamma2 << "," << "," << alpha << ","<< capOmega << "," << beta << std::endl; 

cx_mat lind=com+di1+di2;

return lind;
}


cx_mat Rk4(double g1, double g2, double F, int N, double dt, cx_mat rho)
{
  cx_mat L(double,double,double, int, cx_mat);
  cx_mat K1=dt*L(g1,g2,F,N,rho);
  cx_mat K2=dt*L(g1,g2,F,N,rho+0.5*K1);
  cx_mat K3=dt*L(g1,g2,F,N,rho+0.5*K2);
  cx_mat K4=dt*L(g1,g2,F,N,rho+K3);
  rho=rho+((1.0/6.0)*(K1+(2.0*K2)+(2.0*K3)+K4));
  return rho;

}


int main()
{
ofstream max_energy_out("max_energy_time_g1_0.1_g2_0.1_F_diff.dat");
ofstream max_ergo_out("max_ergotropy_time_g1_0.1_g2_0.1_F_diff.dat");

int N=10; double g1=0.1, g2=0.1, t=150; double w0=1.06;
// Array of F values
vector<double> F_values = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5};


// below for max energy
for (double F : F_values) {

mat I2 = eye(2,2);
mat pauliz={{1,0},{0,-1}}; mat Hb=w0*0.5*(pauliz+I2);

//----eigenvalues of hamiltonian Hb-----------
vec e_n; mat ket_e_n; // eigenvals and vects of H_B
eig_sym(e_n,ket_e_n, Hb);


double t0=0.0;
double dt=0.01;
double l1=(t-t0)/dt;
int m1=round(l1);

// Initial state---------------------
cx_mat rhot=rh0(N);
  
cx_mat rhob0=TrX(rhot, {1}, {10, 2});
vec e_n2; cx_mat ket_e_n2; // eigenvals and vects of H_B
eig_sym(e_n2,ket_e_n2, rhob0);
//------------------------------------
// Init ergotropy and energy----------
//complex<double> E_p = energy_passive(rhob,e_n);

double r0=e_n(0), r1=e_n(1);
double e0=e_n2(0), e1=e_n2(1);
double E_p0=r0*e1+r1*e0;

complex<double> energy0=trace(rhob0*Hb);
complex<double> Ergotropy0=energy0-E_p0;

// ---------------------------------------
//------------solving rk4-----------------------------------------
// Variables to store maximum energy 
double max_energy = real(energy0);

double max_energy_time = t0;

for(int k=1; k<=m1; k++)
{
  t0=t0+dt;  

  rhot=Rk4(g1,g2,F,N,dt,rhot);
    
  cx_mat rhob=TrX(rhot, {1}, {10, 2});


  vec e_n2; cx_mat ket_e_n2; // eigenvals and vects of H_B
  eig_sym(e_n2,ket_e_n2, rhob);
  //------------------------------------
  // Init ergotropy and energy----------
  //complex<double> E_p = energy_passive(rhob,e_n);

  double s0=e_n(0), s1=e_n(1);
  double f0=e_n2(0), f1=e_n2(1);
  double E_p=s0*f1+s1*f0;

  complex<double> energy=trace(rhob*Hb);
  complex<double> Ergotropy=energy-E_p;

      // Update maximum values and corresponding times if new maxima found
      if (real(energy) > max_energy) {
          max_energy = real(energy);
          max_energy_time = t0;
      }

      // Break the loop if energy and ergotropy drop after reaching their maximum values
    if (real(energy) < max_energy) {
        break;
}
}
cout << "For F = " << F << ":" << endl;
cout << "Maximum energy: " << max_energy << " at time: " << max_energy_time << endl;
max_energy_out<< F << '\t' << max_energy << '\t' << max_energy_time << endl;

}

// below for max ergotropy
for (double F : F_values) {

mat I2 = eye(2,2);
mat pauliz={{1,0},{0,-1}}; mat Hb=w0*0.5*(pauliz+I2);

//----eigenvalues of hamiltonian Hb-----------
vec e_n; mat ket_e_n; // eigenvals and vects of H_B
eig_sym(e_n,ket_e_n, Hb);


double t0=0.0;
double dt=0.01;
double l1=(t-t0)/dt;
int m1=round(l1);

// Initial state---------------------
cx_mat rhot=rh0(N);
  
cx_mat rhob0=TrX(rhot, {1}, {10, 2});
vec e_n2; cx_mat ket_e_n2; // eigenvals and vects of H_B
eig_sym(e_n2,ket_e_n2, rhob0);
//------------------------------------
// Init ergotropy and energy----------
//complex<double> E_p = energy_passive(rhob,e_n);

double r0=e_n(0), r1=e_n(1);
double e0=e_n2(0), e1=e_n2(1);
double E_p0=r0*e1+r1*e0;

complex<double> energy0=trace(rhob0*Hb);
complex<double> Ergotropy0=energy0-E_p0;

// ---------------------------------------
//------------solving rk4-----------------------------------------
// Variables to store maximum energy 
double max_ergotropy = real(Ergotropy0);

double max_ergotropy_time = t0;


for(int k=1; k<=m1; k++)
{
  t0=t0+dt;  

  rhot=Rk4(g1,g2,F,N,dt,rhot);
    
  cx_mat rhob=TrX(rhot, {1}, {10, 2});


  vec e_n2; cx_mat ket_e_n2; // eigenvals and vects of H_B
  eig_sym(e_n2,ket_e_n2, rhob);
  //------------------------------------
  // Init ergotropy and energy----------
  //complex<double> E_p = energy_passive(rhob,e_n);

  double s0=e_n(0), s1=e_n(1);
  double f0=e_n2(0), f1=e_n2(1);
  double E_p=s0*f1+s1*f0;

  complex<double> energy=trace(rhob*Hb);
  complex<double> Ergotropy=energy-E_p;

      // Update maximum values and corresponding times if new maxima found
      if (real(Ergotropy) > max_ergotropy) {
          max_ergotropy = real(Ergotropy);
          max_ergotropy_time = t0;
      }

      // Break the loop if energy and ergotropy drop after reaching their maximum values
    if (real(Ergotropy) < max_ergotropy) {
        break;
}
}
cout << "For F = " << F << ":" << endl;
cout << "Maximum ergotropy: " << max_ergotropy << " at time: " << max_ergotropy_time << endl;
max_ergo_out<< F << '\t' << max_ergotropy << '\t' << max_ergotropy_time << endl;

}

return 0;
}
