#include<iostream>
#include<string>
#include<cmath>
#include<math.h>
#include<complex>

using namespace std;

const float VEV = 246;
const float m_h = 125.2;
const float Gamma_h = 4.07e-3;

const float M_W = 80.3692;
const float M_Z = 91.188;

const float m_e = 0.511e-3;
const float m_muon = 105.66e-3;
const float m_tau = 1.77693;
const float m_u = 2.16e-3;
const float m_d = 4.7e-3;
const float m_c = 1.273;
const float m_s = 93.5e-3;
const float m_t = 172.57;
const float m_b = 4.183;

float M_V[2] = {M_W, M_Z};
float M_f[9] = {m_e,m_muon,m_tau,m_u,m_d,m_c,m_s,m_t,m_b};
int n_cs[9] = {1,1,1,3,3,3,3,3,3};
int g_V = 3;
int g_f = 2;

struct particle {float mass;
        string type;
        int dof;
        int n_cs;
        };

particle higgs, Wboson, Zboson, electron, muon, tau, up, down, charm, strange, top, bottom;

higgs.mass = m_h;
higgs.type = "scalar";
higgs.dof = 3;
higgs.n_cs = 1;

const float C = 0.349;

const float gstar = 106.75;
const float MP = 2.4e18;
const float T0 = 2.725*8.617333262e-5*1e-9;
const float gstars0 = 3.91;
const float gstar0 = 3.38;
const float s0 = 2*pow(M_PI,2)/45*gstars0*pow(T0,3);
const float gS = 1;

string particles[12] = {"h","W","Z","e","mu","tau","u","d","c","s","t","b"};

//-------------------------------------------------------------------

double sigma(double w,double MS, double lambda_HS, string part){
    if (part == "h"){
        float M = m_h;
    }
    return 0.;
}

//-------------------------------------------------------------------
int main(){
    complex<double> c1( 4.0 , 5.0 );
    cout << "Specifying initial real & imaginary parts,"
        << "c1 = " << c1 << endl;
    cout<<c1.real()<<endl;

    cout<<higgs.mass<<endl;
    return 0;
}