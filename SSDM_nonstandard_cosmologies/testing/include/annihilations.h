#pragma once

#include <string>
#include <functional>

#include "constants.h"

class sigmaA{

    public:
        long double MS;
        long double lHS;
        particle part;
        sigmaA() = default;
        sigmaA(long double MS,long double lHS,particle part);      
        std::function<long double(long double)> sigmaXX;

};

class sigmaS: public sigmaA{
    public:
        sigmaS() = default;
        sigmaS(long double MS,long double lHS,particle part);
        std::function<long double(long double)> sigmaSS;
        std::function<long double(long double,bool)> sigmav;
};

