#pragma once

#include <iostream>
#include <string>
#include <vector>

using namespace std;
/************************************************************************/
/* Useful fuctions: linear interpolation, finite-difference derivatives */
/************************************************************************/

//Linear interpolation function
long double interp(long double x, const vector<long double> &xData,
              const vector<long double> &yData, bool extrapolate);

//Array of slopes using Finite-difference formulas
vector<long double> slopearray(const vector<long double>& xData,
                          const vector<long double>& yData);
/***********************************************************/
/* g*(S): Effective number of degrees of freedom in the SM */
/***********************************************************/
//Read gstar(S) data from various .tab files in the gstar folder, and compute
//their finite-difference derivatives as a function of T
extern vector<long double> Tvec;

extern vector<long double> gstarvec;
extern vector<long double> dlngstardlnTvec;

extern vector<long double> gstarSvec;
extern vector<long double> dlngstarSdlnTvec;

void Read_gstar(const string& choice, const string& gstarpath);
//g*
long double gstar(long double T);

//g*S
long double gstarS(long double T);

//dlng*S/dlnT
long double dlngstarSdlnT(long double T);

//dlng*/dlnT
long double dlngstardlnT(long double T);

/***************************************************************************/
/* Energy density, comoving entropy, and Hubble rate in the Visible sector */
/***************************************************************************/

//Rho Visible
long double RhoVisible(long double T);

//Comoving entropy of the Visible sector
long double EntropyVisible(long double T);

//Hubble rate
long double Hubble(long double T);

//(H / Hbar) to account for varying gstarS only in the Visible sector
long double HoverHbarVisible(long double T);
