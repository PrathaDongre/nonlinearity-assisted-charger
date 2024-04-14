

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

/*complex<double> energy_passive(cx_mat rhob, cx_vec e_n)
{
 cx_vec r_n; cx_mat ket_r_n; // eigenvalues and eigenvectors of rho
 eig_gen(r_n,ket_r_n, rhob);
 complex<double> E_passive = r_n(1)*e_n(0)+r_n(0)*e_n(1);
 return E_passive;
}*/

cx_mat L(double g1, double g2, double F, int N, cx_mat rho0){
  double w0=1.0;
  double w1=1.0;
  mat I2 = eye(2,2);

  mat paulix={{0,1},{1,0}}; cx_mat pauliy={{0,-im},{im,0}}; mat pauliz={{1,0},{0,-1}};
  cx_mat splus=(paulix+im*pauliy)/2.0; cx_mat sminus=(paulix-im*pauliy)/2.0;
  
  vec eigenval_Ha; 
  mat eigenvec_Ha;
  mat a=aa(N); mat adag=aadagger(N);

  // --------hamiltonian of Anharmonic Oscillator

  mat Ha = w0*a*adag - w1*(a*adag)*(a*adag);

  eig_sym(eigenval_Ha,eigenvec_Ha,Ha);

  //------------------- commutator-----------------------------------------------------------

  cx_mat f1=g1*(tensor(a,splus)+tensor(adag,sminus))+g2*(tensor(a*a,splus)+tensor(adag*adag,sminus))+F*(tensor(a,I2)+tensor(adag,I2));
  cx_mat com=-im*(f1*rho0-rho0*f1);

  //-----------------dissipator--------------------------------------------------------
  //double nb=1/(exp(w0/T)-1.0);

  double nb=0.0,gamma=w0;

  // ------------first term of dissipators---------------------------------------------
  int n = 0;
  cx_mat transition_down = sqrt(n+1) * eigenvec_Ha.col(n) * trans(eigenvec_Ha.col(n+1));
  cx_mat transition_up = sqrt(n+1) * eigenvec_Ha.col(n+1) * trans(eigenvec_Ha.col(n));

  cx_mat di1 = zeros<cx_mat>(rho0.n_rows, rho0.n_cols);
  cx_mat di2 = zeros<cx_mat>(rho0.n_rows, rho0.n_cols);

  //-----------------looping over all m for terms of dissipators---------------------------
  for (int n = 0; n <=N-1; n++)
  {

    cx_mat transition_down = sqrt(n+1) * eigenvec_Ha.col(n) * trans(eigenvec_Ha.col(n+1));
    cx_mat transition_up = sqrt(n+1) * eigenvec_Ha.col(n+1) * trans(eigenvec_Ha.col(n));

    cx_mat di1m = tensor(transition_down,I2)*rho0*tensor(transition_up,I2)-0.5*(tensor(transition_up,I2)*tensor(transition_down,I2)*rho0+rho0*tensor(transition_up,I2)*tensor(transition_down,I2));
    cx_mat di2m = tensor(transition_up,I2)*rho0*tensor(transition_down,I2)-0.5*(tensor(transition_down,I2)*tensor(transition_up,I2)*rho0+rho0*tensor(transition_down,I2)*tensor(transition_up,I2));


    di1 += di1m;
    di2 += di2m;

    di1 *= gamma*(nb+1.0);
    di2 *= gamma*(nb);
  }

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
ofstream out("energy_anharmonic_test.dat");
ofstream ergo_out("ergotropy_anharmonic_test.dat");

int N=10; double g1=0.1, g2=0.1, F=0.5, t=150; double w0=1.0;

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

out << t0 << '\t' << real(energy0) << endl; 
ergo_out << t0 << '\t' << real(Ergotropy0) <<endl;
cout << t0 << '\t' << energy0 << endl;

//------------solving rk4-----------------------------------------

for(int k=1; k<=m1; k++)
{
t0=t0+dt;  

rhot=Rk4(g1,g2,F,N,dt,rhot);
  
cx_mat rhob=TrX(rhot, {1}, {10, 2});


vec e_n2; cx_mat ket_e_n2; // eigenvals and vects of H_B
eig_sym(e_n2,ket_e_n2, rhob);


double s0=e_n(0), s1=e_n(1);
double f0=e_n2(0), f1=e_n2(1);
double E_p=s0*f1+s1*f0;

complex<double> energy=trace(rhob*Hb);
complex<double> Ergotropy=energy-E_p;

out << t0 << '\t' << real(energy) << endl; 
ergo_out << t0 << '\t' << real(Ergotropy) <<endl;
cout << t0 << '\t' << energy << endl;
}


return 0;
}
