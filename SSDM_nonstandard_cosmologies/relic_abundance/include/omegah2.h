#pragma once

#include "particle.h"

long double CollisionTerm(long double T,long double MS, long double lHS, Particle part);
long double HiggsDecay(long double T,long double MS,long double lHS);
long double Abundance(long double MS, long double lHS,long double TR);