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

sigmaA::sigmaA(long double MS_, long double lHS_, particle p)
    : MS(MS_), lHS(lHS_), part(p)
{
    sigmaXX = [MS_, lHS_, p](long double s) {
        std::string name = p.name;
        long double M = p.mass;
        long double nc = p.nc;
        int ID = p.ID;
        long double sigma;
        long double A;
        long double B;
        long double C;

        switch(ID){
            case 0:{
                
                long double xi = sqrt((s-4*M*M)*(s-4*MS_*MS_))/(s-2*M*M);
                long double Fxi = atanh(xi)/xi;

                A = lHS_*lHS_/(8*M_PI*s);
                B = sqrt((s-4*M*M)/(s-4*MS_*MS_));
                C = powl((s+2*M*M)/(s-M*M),2)-(16*lHS_*vev*vev/(s-2*M*M))*((s+2*M*M)/(s-M*M))*Fxi+32*lHS_*lHS_*powl(vev,4)/((s-2*M*M)*(s-2*M*M))*(1/(1-xi*xi)+Fxi);
                sigma = A*B*C;
                break;}
            case 1:{
                A = lHS_*lHS_/(4*M_PI*s);
                B = sqrt((s-4*M*M)/(s-4*MS_*MS_));
                C = (s*s-4.0*s*M*M+12.0*M*M*M*M)/((s-MH*MH)*(s-MH*MH)+(MH*Gammah)*(MH*Gammah));
                sigma = A*B*C;
                break;}
            case 2:{
                A = lHS_*lHS_/(8*M_PI*s);
                B = sqrt((s-4*M*M)/(s-4*MS_*MS_));
                C = (s*s-4.0*s*M*M+12.0*M*M*M*M)/((s-MH*MH)*(s-MH*MH)+(MH*Gammah)*(MH*Gammah));
                sigma = A*B*C;
                break;}
            default:{
                A = lHS_*lHS_*nc/(2*M_PI*s);
                B = pow(s-4*M*M,1.5)/sqrt(s-4*MS_*MS_);
                C = M*M/((s-MH*MH)*(s-MH*MH)+(MH*Gammah)*(MH*Gammah));
                sigma = A*B*C;}
        }

        return sigma; 
    };
}

sigmaS::sigmaS(long double MS_, long double lHS_, particle p):sigmaA(MS_,lHS_,p){
    sigmaSS = [MS_,lHS_,p,this](long double s) {
        std::string name = p.name;
        long double M = p.mass;
        long double nc = p.nc;
        int ID = p.ID;
        long double sigma;
        long double k;
        switch(ID){
            case 0:{
                k = (s-4*MS_*MS_)/(s-4*M*M);
                sigma = k*this->sigmaXX(s);
                break;}
            case 1:{
                k = 0.5*(s-4*MS*MS)/(9.0L*(s-4*M*M));
                sigma = k*this->sigmaXX(s);
                break;}
            case 2:{
                k = (s-4*MS*MS)/(9.0L*(s-4*M*M));
                sigma = k*this->sigmaXX(s);
                break;}
            default:{
                k = 0.5*(s-4*MS*MS)/((s-4*M*M)*4.0L*nc*nc);
                sigma = k*this->sigmaXX(s);}
        };
        return sigma;
    };
    sigmav = [MS_,lHS_,p,this](long double T, bool getIntegrand){
        long double sv;
        // long double gi = 1.0L;  // g.d.l de la partícula de DM, se pueden cambiar según el modelo
        // auto integrand = [MS_,lHS_,p,T,this](long double s){
        //    return 2*this->sigmaSS(s)*sqrt(s)*(s-4*MS_*MS_)*boost::math::cyl_bessel_k(1,sqrt(s)/T);
        // };
        long double m = max(MS_,p.mass);
        auto integrand = [MS_,lHS_,p,T,this](long double s){
            // Chequeos físicos
            if (s <= 4*MS_*MS_ || s <= 4*p.mass*p.mass) return 0.0L;
            // Chequeo para el caso 0
            if (p.ID == 0) {
                long double num = (s-4*p.mass*p.mass)*(s-4*MS_*MS_);
                long double den = s-2*p.mass*p.mass;
                if (num < 0 || fabs(den) < 1e-10) return 0.0L;
                long double xi = sqrt(num)/den;
                if (fabs(xi) >= 1.0) return 0.0L;
            }
            // long double Bess = boost::math::cyl_bessel_k(1,sqrt(s)/T)*1/(sqrt(s))/pow(boost::math::cyl_bessel_k(2,p.mass/T),2);
            long double Bess = sqrt(2*p.mass*p.mass/(M_PI*T*sqrt(s)))*exp(-(sqrt(s)/T-2*p.mass/T))*K1pol(T/sqrt(s))/pow(K2pol(T/p.mass),2);
            long double val = this->sigmaSS(s)*(s-4*p.mass*p.mass)*sqrt(s)/(8*pow(p.mass,4)*T)*Bess;
            if (std::isnan(val) || std::isinf(val)) return 0.0L;
            return val;
        };

        auto integrand_w = [=](long double w) -> long double {
            if(w==M_PI/4 || w==M_PI/2) return 0.;
            return integrand(4.0L*m*m*powl(tanl(w),2))*4*m*m/powl(cosl(w),2);
        };

        long double s60 = 4.0*(60.0*60.0+m*m);
        long double s65 = 4.0*(65.0*65.0+m*m);
        
        auto integrand_u = [&](double u) -> double {
            // You may need to define IntegrandUParams and params_void appropriately
            // For now, this is a placeholder for your logic
            // IntegrandUParams* params = static_cast<IntegrandUParams*>(params_void);
            if(u==0. || u==1.) return 0.;
            double z = 1 - u*u;
            double MM = 2*m;
            long double s = (MM-3*T*log(z))*(MM-3*T*log(z));
            return integrand(s)*2*(MM-3*T*log(z))*6*T*u/z;
        };

        sv = simpson(integrand_u,0.,1.,1e-3,nullptr);
        // if ((p.ID > 2 && p.ID < 11) && (MS<=MH/2)){
        //     // double error;
        //     // long double sv1 = tanh_sinh<long double>().integrate(integrand,4*m*m+1e-15,s60);
        //     // long double sv2 = gauss_kronrod<double, 60>::integrate(integrand,s60,s65,5,1e-14,&error);//tanh_sinh<long double>().integrate(integrand,s60,s65);
        //     // long double sv3 = exp_sinh<long double>().integrate(integrand,s65,INFINITY);
        //     // sv = sv1 + sv2 + sv3;

        //     sv = simpson(integrand_u,0.,1.,1e-3,nullptr);
        // }

        // else {sv = exp_sinh<long double>().integrate(integrand,4*m*m+1e-15,INFINITY);};

        if (getIntegrand){
            ofstream file;
            string filename;
            filename = "integrand_mine_" + p.name + ".txt";

            file.open(filename);
            if (!file.is_open()) 
                throw std::runtime_error("No se pudo abrir el archivo: "+filename);

            file << "# Integrand" << endl;
            file << "# lambda_HS = "<<lHS<<endl;
            file << "# Incoming particle = "<<p.name<<endl;
            file << "# ================================================="<<endl;
            file << "# s [GeV²]\t\tI(s) [GeV^-4]"<<endl;

            long double s60 = 4*(p.mass*p.mass+60.0*60.0);
            long double s65 = 4*(p.mass*p.mass+65.0*65.0);
            long double s_;
            int N = 1001;
            long double delta = (s65 - s60)/(N-1);
            for (int i=0;i<N;i++){
                s_ = s60 + i*delta;
                file << s_ << "\t" << integrand(s_) << endl;
            }

        }

        // sv = simpson(integrand_u, 0., 1.,eps,NULL);

        // sv = tanh_sinh<long double>().integrate(integrand_w,M_PI/4+1e-10,M_PI/2-1e-10);

        return sv;
    };
}
