#pragma once

#include "particle.h"

long double crossSection(long double s, long double MS, long double lambdaHS, Particle part);
long double decayRate(long double MS, long double lambdaHS);
long double neq(long double T,Particle part);
long double sigmav(long double T,long double MS,long double lHS,Particle part);
long double sigmavT(long double T,long double MS,long double lHS);