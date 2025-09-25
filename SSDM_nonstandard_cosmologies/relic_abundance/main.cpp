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
#include "DMproduction.h"
#include "functions.h"
#include "constants.h"
#include "omegah2.h"
#include "particle.h"
#include "annihilations.h"
#include "prod.h"


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
    Read_gstar("standard","/home/santiago/FreezeIn");

//     long double MS = pow(10,2.530612);
//     int mass_name = MS;
//     long double lHS = pow(10,-11.038680);
//     long double T = 1e6;

//     vector<Particle> channels;
//     channels.push_back(Particle("H",MH,1,1));
//     channels.push_back(Particle("W",MW,3,1));
//     channels.push_back(Particle("Z",MZ,3,1));
//     channels.push_back(Particle("t",Mt,2,3));
//     channels.push_back(Particle("b",Mb,2,3));

//     for(auto &p:channels){
//         string filename = "sigmav"+p.name+to_string(mass_name)+".csv";
//         ofstream File(filename);
//         File<<"Dark Matter mass = "<<MS<<" GeV\n";
//         File<<"Coupling = "<<lHS;
//         File<<"T in GeV, sigmav in GeV^-2\n";
//         File<<"Channel: "<<p.name<<"\n";
//         File<<"#############################\n";

//         File<<"T,sigmav\n";
//         int N = 1000;
//         long double i1 = log10l(T0);
//         long double i2 = log10l(T);

//         for(int i = 0;i < N;i++){
//         long double n = i1 + i*(i2-i1)/(N-1);
//         double temp = pow(10.,n);
//         double sv = sigmav(temp,MS,lHS,p);
//         File<<temp<<","<<sv<<"\n";
// }
//     }

//     ofstream File("sigmavTotal.csv");
//     File<<"Dark Matter mass = "<<MS<<" GeV\n";
//     File<<"Coupling = "<<lHS;
//     File<<"T in GeV, sigmav in GeV^-2\n";
//     File<<"Channel: all\n";
//     File<<"#############################\n";

//     File<<"T,sigmav\n";
//     int N = 1000;
//     long double i1 = log10l(T0);
//     long double i2 = log10l(T);

//     for(int i = 0;i < N;i++){
//     long double n = i1 + i*(i2-i1)/(N-1);
//     double temp = pow(10.,n);
//     double sv = sigmavT(temp,MS,lHS);
//     File<<temp<<","<<sv<<"\n";}
//     return 0;

// vector <long double> Ms,Gamma;
// long double lHS = 1e-11;

// long double dM = log10(MH/2)/49;

// for(int i=0;i<50;i++){
//     Ms.push_back(pow(10,i*dM));
// }

// for(int i=0;i<50;i++){
//     Gamma.push_back(decayRate(Ms[i],lHS));
// }

// ofstream File("Gamma_mine.csv");
// File<<"M,Gamma\n";
// for(int i=0;i<50;i++){
//    File<<Ms[i]<<","<<Gamma[i]<<"\n"; 
// }

//Particle part = Particle("b",Mb,2,3);


long double p = atof(argv[1]);
long double MS = atof(argv[2]);
long double s = 4*(p*p+MS*MS);
long double lHS = atof(argv[3]);

cout<<p<<"\t"<<MS<<"\t"<<lHS<<endl;
cout<<"---------------------------------"<<endl;
cout<<"x1,x1->Z,Z"<<"\t"<<sigmaZZ(s,MS,lHS)/pb<<" pb"<<endl;
cout<<"x1,x1->W+,W-"<<"\t"<<sigmaWW(s,MS,lHS)/pb<<" pb"<<endl;
cout<<"x1,x1->b,B"<<"\t"<<sigmaff(s,MS,lHS,Mb,3)/pb<<" pb"<<endl;
cout<<"x1,x1->t,T"<<"\t"<<sigmaff(s,MS,lHS,Mt,3)/pb<<" pb"<<endl;
cout<<"x1,x1->c,C"<<"\t"<<sigmaff(s,MS,lHS,Mc,3)/pb<<" pb"<<endl;
cout<<"x1,x1->s,S"<<"\t"<<sigmaff(s,MS,lHS,Ms,3)/pb<<" pb"<<endl;
cout<<"x1,x1->u,U"<<"\t"<<sigmaff(s,MS,lHS,Mu,3)/pb<<" pb"<<endl;
cout<<"x1,x1->d,D"<<"\t"<<sigmaff(s,MS,lHS,Md,3)/pb<<" pb"<<endl;
cout<<"x1,x1->e,E"<<"\t"<<sigmaff(s,MS,lHS,Me,1)/pb<<" pb"<<endl;
cout<<"x1,x1->m,M"<<"\t"<<sigmaff(s,MS,lHS,Mmu,1)/pb<<" pb"<<endl;  // muon
cout<<"x1,x1->l,L"<<"\t"<<sigmaff(s,MS,lHS,Mta,1)/pb<<" pb"<<endl;  // tau
cout<<"x1,x1->h,h"<<"\t"<<sigmahh(s,MS,lHS)/pb<<" pb"<<endl;
cout<<"---------------------------------"<<endl;
cout<<"#################################"<<endl;

Particle Z = Particle("Z",MZ,3,1);
Particle W = Particle("W",MW,3,1);
Particle bottom = Particle("b",Mb,2,3);
Particle top = Particle("t",Mt,2,3);
Particle charm = Particle("c",Mc,2,3);
Particle strange = Particle("s",Ms,2,3);
Particle up = Particle("u",Mu,2,3);
Particle down = Particle("u",Md,2,3);
Particle electron = Particle("e",Me,2,1);
Particle muon = Particle("mu",Mmu,2,1);
Particle tau = Particle("tau",Mta,2,1);
Particle higgs = Particle("h",MH,1,1);

cout<<"Z,Z->x1,x1"<<"\t"<<crossSection(s,MS,lHS,Z)/pb<<" pb"<<endl;
cout<<"W+,W-->x1,x1"<<"\t"<<crossSection(s,MS,lHS,W)/pb<<" pb"<<endl;
cout<<"b,B->x1,x1"<<"\t"<<crossSection(s,MS,lHS,bottom)/pb<<" pb"<<endl;
cout<<"t,T->x1,x1"<<"\t"<<crossSection(s,MS,lHS,top)/pb<<" pb"<<endl;
cout<<"c,C->x1,x1"<<"\t"<<crossSection(s,MS,lHS,charm)/pb<<" pb"<<endl;
cout<<"s,S->x1,x1"<<"\t"<<crossSection(s,MS,lHS,strange)/pb<<" pb"<<endl;
cout<<"u,U->x1,x1"<<"\t"<<crossSection(s,MS,lHS,up)/pb<<" pb"<<endl;
cout<<"d,D->x1,x1"<<"\t"<<crossSection(s,MS,lHS,down)/pb<<" pb"<<endl;
cout<<"e,E->x1,x1"<<"\t"<<crossSection(s,MS,lHS,electron)/pb<<" pb"<<endl;
cout<<"m,M->x1,x1"<<"\t"<<crossSection(s,MS,lHS,muon)/pb<<" pb"<<endl;  // muon
cout<<"l,L->x1,x1"<<"\t"<<crossSection(s,MS,lHS,tau)/pb<<" pb"<<endl;  // tau
cout<<"h,h->x1,x1"<<"\t"<<crossSection(s,MS,lHS,higgs)/pb<<" pb"<<endl;
cout<<"---------------------------------"<<endl;
cout<<"#################################"<<endl;

cout<<"Z,Z->x1,x1"<<"\t"<<sigmaZZSS(s,MS,lHS)/pb<<" pb"<<endl;
cout<<"W+,W-->x1,x1"<<"\t"<<sigmaWWSS(s,MS,lHS)/pb<<" pb"<<endl;
cout<<"b,B->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Mb,3)/pb<<" pb"<<endl;
cout<<"t,T->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Mt,3)/pb<<" pb"<<endl;
cout<<"c,C->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Mc,3)/pb<<" pb"<<endl;
cout<<"s,S->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Ms,3)/pb<<" pb"<<endl;
cout<<"u,U->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Mu,3)/pb<<" pb"<<endl;
cout<<"d,D->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Md,3)/pb<<" pb"<<endl;
cout<<"e,E->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Me,1)/pb<<" pb"<<endl;
cout<<"m,M->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Mmu,1)/pb<<" pb"<<endl;  // muon
cout<<"l,L->x1,x1"<<"\t"<<sigmaffSS(s,MS,lHS,Mta,1)/pb<<" pb"<<endl;  // tau
cout<<"h,h->x1,x1"<<"\t"<<sigmahhSS(s,MS,lHS)/pb<<" pb"<<endl;
cout<<"---------------------------------"<<endl;


return 0;

}
