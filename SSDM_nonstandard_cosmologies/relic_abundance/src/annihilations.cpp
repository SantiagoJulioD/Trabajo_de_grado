#include <iostream>
#include <cmath>
#include <string>
#include <complex>
#include <vector>
#include "constants.h"
#include "annihilations.h"
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

long double sigmaZZ(long double s,long double MS,long double lHS){
    long double sigma;
    long double A = lHS*lHS/(8*M_PI*s);
    long double B = sqrt((s-4*MZ*MZ)/(s-4*MS*MS));
    long double C = (s*s-4.0*s*MZ*MZ+12.0*MZ*MZ*MZ*MZ)/((s-MH*MH)*(s-MH*MH)+(MH*Gammah)*(MH*Gammah));
    sigma = A*B*C;
    return sigma;
}

long double sigmaWW(long double s,long double MS,long double lHS){
    long double sigma;
    long double A = lHS*lHS/(4*M_PI*s);
    long double B = sqrt((s-4*MW*MW)/(s-4*MS*MS));
    long double C = (s*s-4.0*s*MW*MW+12.0*MW*MW*MW*MW)/((s-MH*MH)*(s-MH*MH)+(MH*Gammah)*(MH*Gammah));
    sigma = A*B*C;
    return sigma;
}

long double sigmaff(long double s,long double MS,long double lHS,long double mf,long double nc){
    long double sigma;
    long double A = lHS*lHS*nc/(2*M_PI*s);
    long double B = pow(s-4*mf*mf,1.5)/sqrt(s-4*MS*MS);
    long double C = mf*mf/((s-MH*MH)*(s-MH*MH)+(MH*Gammah)*(MH*Gammah));
    sigma = A*B*C;
    return sigma;
}

long double sigmahh(long double s,long double MS,long double lHS){
    long double sigma;
    long double xi = sqrt((s-4*MH*MH)*(s-4*MS*MS))/(s-2*MH*MH);
    long double Fxi = atanh(xi)/xi;

    long double A = lHS*lHS/(8*M_PI*s);
    long double B = sqrt((s-4*MH*MH)/(s-4*MS*MS));
    long double C = powl((s+2*MH*MH)/(s-MH*MH),2)-(16*lHS*vev*vev/(s-2*MH*MH))*((s+2*MH*MH)/(s-MH*MH))*Fxi+32*lHS*lHS*powl(vev,4)/((s-2*MH*MH)*(s-2*MH*MH))*(1/(1-xi*xi)+Fxi);
    sigma = A*B*C;
    return sigma;
}