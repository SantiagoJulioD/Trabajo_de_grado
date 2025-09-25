#pragma once

#include <iostream>
#include <string>
#include <cmath>

using namespace std;

class Particle{
    private: 
        
    public:
        Particle(string,long double,int,int);
        string name;
        long double mass;
        int dof;
        int nc;

};