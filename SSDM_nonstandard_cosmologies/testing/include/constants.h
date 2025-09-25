#pragma once

#include <string>

// Constantes físicas globales (tipo long double para mayor precisión)
inline const long double vev         = 246.0L;
inline const long double Gammah      = 4.07e-3L;
inline const long double branchRatio = 1.0L - 5.792e-1L - 6.240e-2L - 2.165e-4L - 2.876e-4L;

inline const long double geffIni = 106.75L;
inline const long double MPl     = 2.435323077e18L;

inline const long double T0       = 2.3482233139345615e-13L;
inline const long double s0mo     = 3.059942942537058e-35L;
inline const long double rhoc0mo  = 1.1151984914745777e-43L;
inline const long double s0       = 2.2208139221039978e-38L;
inline const long double rhoc0    = 3.677771412179473e-55L;

inline const long double Me = 0.51099895e-3L;
inline const long double Mmu = 105.6583755e-3L;
inline const long double Mta = 1.77686L;
inline const long double Mu = 2.16e-3L;
inline const long double Mb = 4.18L; // 3.172024L; // 
inline const long double Mt = 172.69L;
inline const long double Md = 4.67e-3L;
inline const long double Ms = 93.4e-3L;
inline const long double Mc = 1.27L; // 0.6977778L; // 
inline const long double MZ = 91.1876L;
inline const long double MW = 80.379L;
inline const long double MH = 125.2L;

inline const long double pb = 2.56819e-9L; // GeV^-2

// SM particles

struct particle {
    std::string name;
    long double mass;
    long double dof;
    long double nc;
    int ID;
};

// Higgs boson
inline const particle higgs = {"h", MH, 1.0L, 1.0L, 0};
// Gauge bosons
inline const particle Wboson = {"W", MW, 3.0L, 1.0L, 1};
inline const particle Zboson = {"Z", MZ, 3.0L, 1.0L, 2};
// Leptons
inline const particle electron = {"e", Me, 2.0L, 1.0L, 3};
inline const particle muon = {"m", Mmu, 2.0L, 1.0L, 4};
inline const particle tau = {"l", Mta, 2.0L, 1.0L, 5};
// Quarks
inline const particle up     = {"u", Mu, 2.0L, 3.0L, 6};
inline const particle down   = {"d", Md, 2.0L, 3.0L, 7};
inline const particle charm  = {"c", Mc, 2.0L, 3.0L, 8};
inline const particle strange = {"s", Ms, 2.0L, 3.0L, 9};
inline const particle bottom = {"b", Mb, 2.0L, 3.0L, 10};
inline const particle top    = {"t", Mt, 2.0L, 3.0L, 11};
