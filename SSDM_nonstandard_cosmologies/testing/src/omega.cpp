#include <iostream>
#include <cmath>
#include <string>
#include <complex>
#include <fstream>
#include <vector>
#include "constants.h"
#include "annihilations.h"
#include "micromegas.h"
#include "simpson.h"
#include "omega.h"

//Boost C++ library
#include <boost/math/special_functions/bessel.hpp>/*provides bessel-K function*/
#include <boost/math/quadrature/gauss.hpp>/*provides Gauss-Legendre quadrature*/
#include <boost/math/quadrature/gauss_kronrod.hpp>/*provides Gauss-Legendre quadrature*/
#include <boost/math/quadrature/exp_sinh.hpp>/*provides exp_sinh quadrature*/
#include <boost/math/quadrature/tanh_sinh.hpp>/*provides exp_sinh quadrature*/
#include <boost/math/tools/roots.hpp> /*roots of a function*/
#include <boost/multiprecision/cpp_bin_float.hpp>

//Namespaces
using namespace std;
using namespace boost::math;
using namespace boost::math::quadrature;
using namespace boost::math::tools;
using boost::multiprecision::cpp_bin_float_quad;

long double omegah2(long double MS, long double lHS, long double TR){
    sigmaS channels[12] = {
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



}