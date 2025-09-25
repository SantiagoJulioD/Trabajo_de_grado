/********************************/
/* Standard and Boost Libraries */
/********************************/

//Standard libraries
#include <iostream>/*standard i/o library*/
#include <fstream>/*to read/write file*/
#include <cstdio>/*provides printf*/
#include <string>/*to use string data type*/
#include <cmath>/*provides pow, sqrt, ...*/
#include <vector>/*provides std::vector*/
#include <functional> /*std::function*/
#include <limits> /* std::numeric_limits */
#include <iomanip>
#include "annihilations.h"
#include "prod.h"
#include "constants.h"
#include "simpson.h"


//Boost C++ library
#include <boost/math/special_functions/bessel.hpp>/*provides bessel-K function*/
#include <boost/math/quadrature/gauss.hpp>/*provides Gauss-Legendre quadrature*/
#include <boost/math/quadrature/exp_sinh.hpp>/*provides exp_sinh quadrature*/
#include <boost/math/tools/roots.hpp> /*roots of a function*/

//Namespaces
using namespace std;
using namespace boost::math;
using namespace boost::math::quadrature;
using namespace boost::math::tools;

int main(int argc, char** argv){

long double p = atof(argv[1]);
long double MS = atof(argv[2]);
long double s = 4*(p*p+MS*MS);
long double lHS = atof(argv[3]);

sigmaA ann[12] = {
    sigmaA(MS, lHS, higgs),
    sigmaA(MS, lHS, Wboson),
    sigmaA(MS, lHS, Zboson),
    sigmaA(MS, lHS, electron),
    sigmaA(MS, lHS, muon),
    sigmaA(MS, lHS, tau),
    sigmaA(MS, lHS, up),
    sigmaA(MS, lHS, down),
    sigmaA(MS, lHS, charm),
    sigmaA(MS, lHS, strange),
    sigmaA(MS, lHS, bottom),
    sigmaA(MS, lHS, top)
};

sigmaS prod[12] = {
    sigmaS(MS, lHS, higgs),
    sigmaS(MS, lHS, Wboson),
    sigmaS(MS, lHS, Zboson),
    sigmaS(MS, lHS, electron),
    sigmaS(MS, lHS, muon),
    sigmaS(MS, lHS, tau),
    sigmaS(MS, lHS, up),
    sigmaS(MS, lHS, down),
    sigmaS(MS, lHS, charm),
    sigmaS(MS, lHS, strange),
    sigmaS(MS, lHS, bottom),
    sigmaS(MS, lHS, top)
};

std::vector<long double> Temp;
int npoints = 1000;
long double TR = 1e6;
long double expTR = log10l(TR);
long double expT0 = log10l(T0);
long double del = (expTR-expT0)/(npoints-1);
// long double asym;
Temp.resize(npoints);
for (int i = 0;i < npoints;i++){
    Temp[i] = powl(10,expTR-i*del);
    // asym = M_PI/2*Temp[i]/MS*expl(-2*MS/Temp[i]);
    // cout<<"T = "<<Temp[i]<<", T*K_2(MS/T)² = "<<Temp[i]*pow(boost::math::cyl_bessel_k(2,MS/Temp[i]),2)<<", asym = "<<asym<<endl;
    // cout<<i<<": T = "<<Temp[i]<<", <σv> = "<<prod[0].sigmav(Temp[i])/pb<<endl;
}

cout<<p<<"\t"<<MS<<"\t"<<lHS<<endl;
cout<<"---------------------------------"<<endl;
for (auto &par:ann){
    cout<<"x1,x1->"<<par.part.name<<",¬"<<par.part.name<<"\t"<<par.sigmaXX(s)/pb<<" pb"<<endl;}
cout<<"---------------------------------"<<endl;
cout<<"#################################"<<endl;
for (auto &par:prod){
    cout<<par.part.name<<",¬"<<par.part.name<<"->x1,x1"<<"\t"<<par.sigmaSS(s)/pb<<" pb"<<endl;}
cout<<"---------------------------------"<<endl;

// ofstream file;

// file.open("sigmahh.txt");
// if (!file.is_open()) 
//     throw std::runtime_error("No se pudo abrir el archivo: sigmav.txt");

// //file << std::fixed << std::setprecision(6);

// file << "# Cross section for hh->~x1~x1 process [pb]"<<endl;
// file << "# M_S = "<<MS<<" GeV"<<endl;
// file << "# lambda_HS = "<<lHS<<endl;
// file << "# ================================================="<<endl;
// file << "# s [GeV²]\t\t sigma [pb]"<<endl;

// long double pini = 1000;
// long double pfin = 100000;
// long double deltap = (pfin-pini)/100;

// long double p_,s_;

// for (int i = 0;i < 101;i++){
//     p_ = pini+i*deltap;
//     s_ = 4*(p_*p_+prod[0].part.mass*prod[0].part.mass);
//     file << p_ << "\t"<<prod[0].sigmaSS(s_)/pb<< endl;
// }

ofstream file,filep;

string filename = "sigmav" + to_string((int)MS) + ".txt";

file.open(filename);
if (!file.is_open()) 
    throw std::runtime_error("No se pudo abrir el archivo: sigmav.txt");

//file << std::fixed << std::setprecision(6);

file << "# Thermally-averaged cross section [pb]"<<endl;
file << "# M_S = "<<MS<<" GeV"<<endl;
file << "# lambda_HS = "<<lHS<<endl;
file << "# ================================================="<<endl;
file << "# T [GeV]\t\t";
for (auto &p:prod){
    file << p.part.name << "\t\t";
}
file << endl;
for (int i = 0;i < npoints;i++){
    file << Temp[i] << "\t";
    for (int j=0;j<12;j++){
        file << prod[j].sigmav(Temp[i],0)/pb<<"\t";
    }
    file << endl;
}

filep.open("sigma.txt");
if (!filep.is_open()) 
    throw std::runtime_error("No se pudo abrir el archivo: sigma.txt");

//file << std::fixed << std::setprecision(6);

filep << "# Cross section [pb]"<<endl;
filep << "# M_S = "<<MS<<" GeV"<<endl;
filep << "# lambda_HS = "<<lHS<<endl;
filep << "# ================================================="<<endl;
filep << "# p [GeV]\t\t";
for (auto &p:prod){
    filep << p.part.name << "\t\t";
}
filep << endl;

long double Pmin = 50.0L;
long double Pmax = 100.0L;
int N = 1001;
long double dP = (Pmax-Pmin)/(N-1);
long double PP, ss;

for (int i = 0;i < N;i++){
    PP = Pmin + i*dP;
    filep << PP << "\t";
    for (int j=0;j<12;j++){
        ss = 4*(PP*PP+prod[j].part.mass*prod[j].part.mass);
        filep << prod[j].sigmaSS(ss)/pb<<"\t";
    }
    filep << endl;
}

long double dummy;

dummy = prod[4].sigmav(62.0,1);

cout<<"T0 = "<<T0<<endl;
cout<<"x0 = "<<MS/T0<<endl;

cout<< "1 pb = "<<1/pb<<" 1/GeV²"<<endl;
cout <<"sigma(P=60) = "<<prod[4].sigmaSS(4.0*(60.0*60.0+Mmu*Mmu)) <<endl;

auto f = [](double x) { return std::sin(x); };
double result = simpson(f, 0.0, M_PI, 1e-6, nullptr);
std::cout << "Integral de sin(x) de 0 a pi: " << result << std::endl;
return 0;

return 0;

}
