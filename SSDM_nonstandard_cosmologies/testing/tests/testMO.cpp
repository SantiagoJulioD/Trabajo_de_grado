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
#include <array>
#include "micromegas.h"
#include "simpson.h"
#include "constants.h"


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

// Definimos la función a integrar
double myFunc(double x) {
    return std::sin(x);
}

int main(){

    int err = 0;
    double result = simpson(myFunc, 0.0, M_PI, 1e-6, &err);
    std::cout << "Integral de sin(x) de 0 a pi: " << result << std::endl;

    vector <double> xs,xss,I0s,I0ss,I1s,I1ss,K0s,K0ss,K1s,K1ss,K2s,K2ss,K1p,K2p;
    int N = 1000;
    double xi = 1e-10;
    double xf = 3.0;
    double xff = 3.0;

    double delta = (xf-xi)/(N-1);
    double deltaExp = (log10(xff)-log10(xi))/(N-1);

    xs.resize(N);
    xss.resize(N);

    I0s.resize(N);
    I0ss.resize(N);
    I1s.resize(N);
    I1ss.resize(N);
    K0s.resize(N);
    K0ss.resize(N);
    K1s.resize(N);
    K1ss.resize(N);
    K2s.resize(N);
    K2ss.resize(N);
    K1p.resize(N);
    K2p.resize(N);
    

    ofstream file1,file2,file3;

    file1.open("bessel.txt");
    if (!file1.is_open()) 
        throw std::runtime_error("No se pudo abrir el archivo: bessel.txt");

    file1 << "# Modified Bessel functions" << endl;
    file1 << "# ===================================" << endl;
    file1 << "# x\t\tI0_ap\t\tI0_ex\t\tI1_ap\t\tI1_ex\t\tK0_ap\t\tK0_ex\t\tK1_ap\t\tK1_ex\t\tK2_ap\t\tK2_ex"<<endl;

    file2.open("bessel10.txt");
    if (!file2.is_open()) 
        throw std::runtime_error("No se pudo abrir el archivo: bessel10.txt");

    file2 << "# Modified Bessel functions" << endl;
    file2 << "# ===================================" << endl;
    file2 << "# x\t\tI0_ap\t\tI0_ex\t\tI1_ap\t\tI1_ex\t\tK0_ap\t\tK0_ex\t\tK1_ap\t\tK1_ex\t\tK2_ap\t\tK2_ex"<<endl;

    file3.open("besselapprox.txt");
    if (!file3.is_open())
        throw std::runtime_error("No se pudo abrir el archivo: besselapprox.txt");
    file3 << "# Modified Bessel functions (approximation)" << endl;
    file3 << "# ===================================" << endl;
    file3 << "# x1\t\tK1\t\tK2\t\tx2\t\tK1\t\tK2" << endl;

    for(int i=0;N-1;i++){
        xs[i] = xi + i*delta;
        xss[i] = pow(10,log10(xi) + i*deltaExp);

        I0s[i] = bessI0(xs[i]);
        I1s[i] = bessI1(xs[i]);
        K0s[i] = bessK0(xs[i]);
        K1s[i] = bessK1(xs[i]);
        K2s[i] = bessK2(xs[i]);

        I0ss[i] = bessI0(xss[i]);
        I1ss[i] = bessI1(xss[i]);
        K0ss[i] = bessK0(xss[i]);
        K1ss[i] = bessK1(xss[i]);
        K2ss[i] = bessK2(xss[i]);
        
        file1 << xs[i] << "\t" 
              << I0s[i] << "\t" 
              << boost::math::cyl_bessel_i(0,xs[i]) << "\t"
              << I1s[i] << "\t" 
              << boost::math::cyl_bessel_i(1,xs[i]) << "\t"
              << K0s[i] << "\t" 
              << boost::math::cyl_bessel_k(0,xs[i]) << "\t"
              << K1s[i] << "\t" 
              << boost::math::cyl_bessel_k(1,xs[i]) << "\t"
              << K2s[i] << "\t" 
              << boost::math::cyl_bessel_k(2,xs[i]) << endl;

        file2 << xss[i] << "\t" 
              << I0ss[i] << "\t" 
              << boost::math::cyl_bessel_i(0,xss[i]) << "\t"
              << I1ss[i] << "\t" 
              << boost::math::cyl_bessel_i(1,xss[i]) << "\t"
              << K0ss[i] << "\t" 
              << boost::math::cyl_bessel_k(0,xss[i]) << "\t"
              << K1ss[i] << "\t" 
              << boost::math::cyl_bessel_k(1,xss[i]) << "\t"
              << K2ss[i] << "\t" 
              << boost::math::cyl_bessel_k(2,xss[i]) << endl;

        file3 << xs[i] << "\t"
              << sqrt(M_PI/2/xs[i])*exp(-xs[i])*K1pol(1/xs[i]) << "\t"
              << sqrt(M_PI/2/xs[i])*exp(-xs[i])*K2pol(1/xs[i]) << "\t"
              << xss[i] << "\t"
              << sqrt(M_PI/2/xss[i])*exp(-xss[i])*K1pol(1/xss[i]) << "\t"
              << sqrt(M_PI/2/xss[i])*exp(-xss[i])*K2pol(1/xss[i]) << endl;
    }

    return 0;
}