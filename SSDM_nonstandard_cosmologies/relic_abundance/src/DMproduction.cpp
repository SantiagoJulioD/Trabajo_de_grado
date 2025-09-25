#include <iostream>
#include <cmath>
#include <string>
#include <complex>
#include <vector>
#include "constants.h"
#include "DMproduction.h"
#include "particle.h"

//Boost C++ library
#include <boost/math/special_functions/bessel.hpp>/*provides bessel-K function*/
#include <boost/math/quadrature/gauss.hpp>/*provides Gauss-Legendre quadrature*/
#include <boost/math/quadrature/exp_sinh.hpp>/*provides exp_sinh quadrature*/
#include <boost/math/quadrature/tanh_sinh.hpp>/*provides exp_sinh quadrature*/
#include <boost/math/tools/roots.hpp> /*roots of a function*/

//Namespaces
using namespace std;
using namespace boost::math;
using namespace boost::math::quadrature;
using namespace boost::math::tools;

long double crossSection(long double s, long double MS, long double lambdaHS, Particle part){
    long double mPart = part.mass;
    long double g = part.dof; 
    long double n = part.nc;
    int giter = part.dof;
    long double sigma = 0.0;

    switch (giter){
        case 1:  // scalar boson (Higgs)
        if ((s-4*MS*MS)/(s-4*mPart*mPart)>0 && (s-4*mPart*mPart)!=0){
            long double A = 1/(g*g)*(2*lambdaHS)*(2*lambdaHS)/(32*M_PI*s);
            long double B = sqrt((s-4*MS*MS)/(s-4*mPart*mPart));
            std::complex<long double> aux (s-mPart*mPart,mPart*Gammah);
            std::complex<long double> C;
            C = 2*lambdaHS+3*mPart*1.0/aux-4*vev*vev*2*lambdaHS/(s-2*mPart*mPart);
            sigma = A*B*norm(C);
        }
        else{
            sigma = 0;
        }
        break;
        case 2:  // fermion
        if ((s-4*MS*MS)*(s-4*mPart*mPart)>0){
            long double A = (2*lambdaHS)*(2*lambdaHS)/(16*M_PI*s*n)/(g*g);
            long double B = sqrt((s-4*MS*MS)*(s-4*mPart*mPart));
            long double C = mPart*mPart/((s-MH*MH)*(s-MH*MH)+(MH*Gammah)*(MH*Gammah));
            sigma = A*B*C;
        }
        else {
            sigma = 0;
        }
        break;
        case 3:   // gauge boson
        if ((s-4*MS*MS)/(s-4*mPart*mPart)>0){
            long double A = 1/(g*g)*(2*lambdaHS)*(2*lambdaHS)/(32*M_PI*s);
            long double B = sqrt((s-4*MS*MS)/(s-4*mPart*mPart));
            long double C = (s*s-4*s*mPart*mPart+12*mPart*mPart*mPart*mPart)/((s-MH*MH)*(s-MH*MH)+(MH*Gammah)*(MH*Gammah));
            sigma = A*B*C;
        }
        else {
            sigma = 0;
        }
        break;
    }
    return sigma;
}

long double decayRate(long double MS, long double lambdaHS){
    long double Gamma;
    if (MS < MH/2){   
        Gamma = (2*lambdaHS)*(2*lambdaHS)*vev*vev/(32*M_PI*MH)*sqrt(1-4*(MS/MH)*(MS/MH));
    }
    else{
        Gamma = 0;
    }
    return Gamma;

}

long double neq(long double T,Particle part){
    long double mass = part.mass;
    long double g = part.dof;

    long double aux;
    long double z = mass/T;
    if(z>3.){
        aux = sqrt(M_PI/(2*z))*exp(-z)*(1+15/(8*z)+15*7/(2*(8*z)*(8*z))); 
    }
    else{
        aux = boost::math::cyl_bessel_k(2, z);
    }

    return T/(2.0*M_PI*M_PI)*g*mass*mass*aux;

}

long double sigmav(long double T,long double MS,long double lHS,Particle part){
    long double mx = part.mass;

    auto integrand_s = [=] (long double s) {
        long double aux;
        long double z = sqrt(s)/T;
        if(z>3.){
            aux = sqrt(M_PI/(2*z))*exp(-z)*(1+3/(8*z)-15/(2*(8*z)*(8*z))); 
        }
        else{
            aux = boost::math::cyl_bessel_k(1, sqrt(s)/T);
        }
            return crossSection(s,MS,lHS,part)*(s-4*mx*mx)*sqrt(s)*aux;
        };

    long double res = exp_sinh<long double>().integrate(integrand_s,4.0L*mx*mx+1e-10,INFINITY);

    return T/(8.0*pow(M_PI,4)*pow(neq(T,part),2))*res;

}

long double sigmavT(long double T,long double MS,long double lHS){
    vector<Particle> SM_particles;
    SM_particles.push_back(Particle("H",MH,1,1));
    SM_particles.push_back(Particle("W",MW,3,1));
    SM_particles.push_back(Particle("Z",MZ,3,1));
    SM_particles.push_back(Particle("e",Me,2,1));
    SM_particles.push_back(Particle("mu",Mmu,2,1));
    SM_particles.push_back(Particle("tau",Mta,2,1));
    SM_particles.push_back(Particle("u",Mu,2,3));
    SM_particles.push_back(Particle("d",Md,2,3));
    SM_particles.push_back(Particle("c",Mc,2,3));
    SM_particles.push_back(Particle("s",Ms,2,3));
    SM_particles.push_back(Particle("t",Mt,2,3));
    SM_particles.push_back(Particle("b",Mb,2,3));

    long double neqTotal = 0;
    long double sv = 0;

    for (auto &p:SM_particles){
        neqTotal += neq(T,p);
        sv += sigmav(T,MS,lHS,p)*neq(T,p)*neq(T,p);
    }
    return sv/(neqTotal*neqTotal);

}