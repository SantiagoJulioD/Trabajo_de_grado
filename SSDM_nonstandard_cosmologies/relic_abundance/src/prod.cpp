#include <iostream>
#include <cmath>
#include <string>
#include <complex>
#include <vector>
#include "constants.h"
#include "annihilations.h"
#include "prod.h"
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

long double sigmaZZSS(long double s,long double MS,long double lHS){
    long double sigma, k;
    k = (s-4*MS*MS)/(s-4*MZ*MZ);
    sigma = sigmaZZ(s,MS,lHS)/9.0L;
    return k*sigma;
}

long double sigmaWWSS(long double s,long double MS,long double lHS){
    long double sigma, k;
    k = 0.5*(s-4*MS*MS)/(s-4*MW*MW);
    sigma = sigmaWW(s,MS,lHS)/9.0L;
    return k*sigma;
}

long double sigmaffSS(long double s,long double MS,long double lHS,long double mf,long double nc){
    long double sigma, k;
    k = 0.5*(s-4*MS*MS)/(s-4*mf*mf);
    sigma = sigmaff(s,MS,lHS,mf,nc)/(4.0L*nc*nc);
    return k*sigma;
}

long double sigmahhSS(long double s,long double MS,long double lHS){
long double sigma, k;
    k = (s-4*MS*MS)/(s-4*MH*MH);
    sigma = sigmahh(s,MS,lHS);
    return k*sigma;
}