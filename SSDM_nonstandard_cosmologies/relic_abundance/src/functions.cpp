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
#include "functions.h"
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

/************************************************************************/
/* Useful fuctions: linear interpolation, finite-difference derivatives */
/************************************************************************/

//Linear interpolation function
long double interp(long double x, const vector<long double> &xData,
              const vector<long double> &yData, bool extrapolate) {

    //Check if x is increasing/decreasing
    bool increasing = xData[1] > xData[0];

    //Perform binary search to find the interval for interpolation
    int mid, low = 0, high = xData.size() - 1;
    while (high - low > 1) {
        
        mid = (low + high)/2;

        if (
                (increasing && x >= xData[mid]) ||
                (!increasing && x <= xData[mid])
           ) {
            low = mid;
        }
        else {
            high = mid;
        }
    }

    //Points on either side (unless beyond ends)
    long double xL = xData[low], yL = yData[low];
    long double xR = xData[low+1], yR = yData[low+1];

    if ( !extrapolate ) {//if beyond ends of array and not extrapolating
       if ( (increasing && x < xL) || (!increasing && x > xL) ) yR = yL;
       if ( (increasing && x > xR) || (!increasing && x < xR) ) yL = yR;
    }

    long double dydx = (yR - yL)/(xR - xL);//gradient
    
    return yL + dydx * ( x - xL );//linear interpolation
}

//Array of slopes using Finite-difference formulas
vector<long double> slopearray(const vector<long double>& xData,
                          const vector<long double>& yData) {

    long double x, y, x1, y1, x2, y2;
    vector<long double> dydxData;

    for (size_t i = 0; i < xData.size(); i++) {

        x = xData[i]; y = yData[i];
        if ( i == 0 ) {
            x1 = xData[i+1]; y1 = yData[i+1]; x2 = xData[i+2]; y2 = yData[i+2];
        }
        else if ( i == xData.size() - 1 ) {
            x1 = xData[i-1]; y1 = yData[i-1]; x2 = xData[i-2]; y2 = yData[i-2];
        }
        else {
            x1 = xData[i-1]; y1 = yData[i-1]; x2 = xData[i+1]; y2 = yData[i+1];
        }

        dydxData.push_back(y*(2.0L*x - (x1 + x2))/(x - x1)/(x - x2) +
                           y1*(x - x2)/(x1 - x)/(x1 - x2) +
                           y2*(x - x1)/(x2 - x)/(x2 - x1));
    }

    return dydxData;
}

/***********************************************************/
/* g*(S): Effective number of degrees of freedom in the SM */
/***********************************************************/

//Define arrays for Temperature T and gstar(S) in the SM

vector<long double> Tvec;

vector<long double> gstarvec;
vector<long double> dlngstardlnTvec;

vector<long double> gstarSvec;
vector<long double> dlngstarSdlnTvec;

//Read gstar(S) data from various .tab files in the gstar folder, and compute
//their finite-difference derivatives as a function of T
void Read_gstar(const string& choice, const string& gstarpath) {
    
    long double r1, r2, r3;
    vector <long double> tempvec;

    //Clear global vectors before filling them
    Tvec.clear(), Tvec.shrink_to_fit();
    gstarvec.clear(), gstarvec.shrink_to_fit();
    gstarSvec.clear(), gstarSvec.shrink_to_fit();
    dlngstardlnTvec.clear(), dlngstardlnTvec.shrink_to_fit();
    dlngstarSdlnTvec.clear(), dlngstarSdlnTvec.shrink_to_fit();

    //Open a file corresponding to the input choice
    string filename;
    if (choice == "standard") { filename = "gstar/std.tab"; }
    else if (choice == "HP_A") { filename = "gstar/HP_A.tab"; }
    else if (choice == "HP_B") { filename = "gstar/HP_B.tab"; }
    else if (choice == "HP_B2") { filename = "gstar/HP_B2.tab"; }
    else if (choice == "HP_B3") { filename = "gstar/HP_B3.tab"; }
    else if (choice == "HP_C") { filename = "gstar/HP_C.tab"; }

    //prepend the path to the gstar folder to the filename
    filename = gstarpath + "/" + filename;

    ifstream file(filename);
    //If a file is open, read line-by-line to extract three values (r1, r2, r3)
    //from each line and store them in the corresponding global vectors
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;//Skip lines with '#'
            istringstream iss(line);
            iss >> r1 >> r2 >> r3;
            Tvec.push_back(r1);
            gstarvec.push_back(r3);
            gstarSvec.push_back(r2);
        }
        file.close();
    }
    else {
        cout << "Unable to open the file:" << filename << endl;
    }

    tempvec = slopearray(Tvec, gstarvec);
    for (size_t i = 0; i < Tvec.size(); i++) {
        dlngstardlnTvec.push_back((Tvec[i]/gstarvec[i])*tempvec[i]);
    }
    tempvec = slopearray(Tvec, gstarSvec);
    for (size_t i = 0; i < Tvec.size(); i++) {
        dlngstarSdlnTvec.push_back((Tvec[i]/gstarSvec[i])*tempvec[i]);
    }

    return;
}
//g*
long double gstar(long double T) {
    return interp(T, Tvec, gstarvec, false);
}

//g*S
long double gstarS(long double T) {
    return interp(T, Tvec, gstarSvec, false);
}

//dlng*S/dlnT
long double dlngstarSdlnT(long double T) {
    return interp(T, Tvec, dlngstarSdlnTvec, false);
}

//dlng*/dlnT
long double dlngstardlnT(long double T) {
    return interp(T, Tvec, dlngstardlnTvec, false);
}

/***************************************************************************/
/* Energy density, comoving entropy, and Hubble rate in the Visible sector */
/***************************************************************************/

//Rho Visible
long double RhoVisible(long double T) {
    return (M_PI*M_PI/30.0L)*gstar(T)*pow(T, 4.0L);
}

//Comoving entropy of the Visible sector
long double EntropyVisible(long double T) {
    return (2.0L*M_PI*M_PI/45.0L)*gstarS(T)*T*T*T;
}

//Hubble rate
long double Hubble(long double T) {
    return sqrt(M_PI*M_PI*gstar(T)/90.0L)*T*T/MPl;
}

//(H / Hbar) to account for varying gstarS only in the Visible sector
long double HoverHbarVisible(long double T) {
    return (1.0L + (1.0L/3.0L)*dlngstarSdlnT(T));
}
