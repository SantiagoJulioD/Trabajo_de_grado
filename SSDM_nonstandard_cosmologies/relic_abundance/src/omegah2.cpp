#include <iostream>/*standard i/o library*/
#include <fstream>/*to read/write file*/
#include <cstdio>/*provides printf*/
#include <string>/*to use string data type*/
#include <cmath>/*provides pow, sqrt, ...*/
#include <vector>/*provides std::vector*/
#include <functional> /*std::function*/
#include <limits> /* std::numeric_limits */
#include "functions.h"
#include "constants.h"
#include "omegah2.h"
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

long double CollisionTerm(long double T,long double MS, long double lHS, Particle part){
    long double mx = part.mass;
    long double gx = part.dof;
    auto integrand_s = [=] (long double s) {
            return crossSection(s,MS,lHS,part)*(s-4*mx*mx)*sqrt(s)*boost::math::cyl_bessel_k(1, sqrt(s)/T);
        };
        
        return (gx*gx)/(32*pow(M_PI,4)) *
               exp_sinh<long double>().integrate(integrand_s,
                                            4.0L*mx*mx+1e-10,
                                            INFINITY);

}

long double HiggsDecay(long double T,long double MS,long double lHS){
    return branchRatio*MH*decayRate(MS,lHS)*T/(2*M_PI*M_PI)*boost::math::cyl_bessel_k(1, MH/T);
}

long double Abundance(long double MS, long double lHS,long double TR){
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

    auto dYdT = [=] (long double T){
        long double term1 = 0.0;
        for(auto &p:SM_particles){
            term1 += CollisionTerm(T,MS,lHS,p);

        }
        return 1/(EntropyVisible(T)*Hubble(T)/HoverHbarVisible(T)*T)*(term1 + HiggsDecay(T,MS,lHS));
    };

    // cout<<dYdT(1e6)<<endl;
    // long double h = 0.67;

    return MS*s0mo*gauss<long double,3001>().integrate(dYdT, T0, TR)/rhoc0mo;

}