
#include <iostream>
#include <string>
#include <cmath>
#include "particle.h"

using namespace std;

Particle::Particle(string text,long double m,int g,int n){
    name = text;
    mass = m;
    dof = g;
    nc = n;
}