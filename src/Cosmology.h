#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <math.h>

inline double betaNormal(double beta, int t, int tSize, double tScale, double(*scaleFactorFunc)(double, int tSize, double tScale))
{
    return beta;
}

inline double betaEff(double beta, int t, int tSize, double tScale, double(*scaleFactorFunc)(double, int tSize, double tScale))
{
    double a;

    a = scaleFactorFunc(t, tSize, tScale) / scaleFactorFunc(0, tSize, tScale);
    return beta * (a*a*a*a*a + a*a*a*a*a*a*a) / 2.0;
}

inline double aMinkowski(double t, int tSize, double tScale)
{
    return 1.0;
}

inline double aRadiation(double t, int tSize, double tScale) // a(t) ~ t^(1/2)
{
    if (t < tSize / 2)
       return sqrt(1.0 + t / tScale);
    else
        return sqrt(1.0 + (tSize - t) / tScale);
}

inline double aMatter(double t, int tSize, double tScale) // a(t) ~ t^(2/3)
{
    if (t < tSize / 2)
        return pow(1 + t / tScale, 2.0 / 3.0);
    else
        return pow(1 + (tSize - t) / tScale, 2.0 / 3.0);
}

inline double aLambda(double t, int tSize, double tScale) // a(t) ~ exp(t)
{
    if (t < tSize / 2)
        return exp(1 + t / tScale);
    else
        return exp(1 + (tSize - t) / tScale);
        
}

#endif
