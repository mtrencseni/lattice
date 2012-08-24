#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>

class Configuration
{
public:
    typedef enum {MINKOWSKI, RADIATION, MATTER, LAMBDA } Cosmology;
    typedef enum {NORMAL, EFFECTIVE} Mode;
    typedef enum {UNITY, RANDOM} Init;

    int         tSize;
    int         sSize;
    Cosmology   cosmology;
    Mode        mode;
    Init        init;
    std::string filename;
    double      betaMin;
    double      betaMax;
    double      dBeta;
    int         iterations;
    double      tScale;
};

#endif
