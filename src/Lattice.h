#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <fstream>
#include "mersenne/Randomc.h"
#include "Configuration.h"
#include "Cosmology.h"

class Lattice
{
public:
    Lattice();
    void            Run(Configuration& config);

private:
    void            MoveUp(int x[], int d);
    void            MoveDown(int x[], int d);
    void            InitRand(double r);
    void            InitLattice();
    void            CalcScaleFactor();
    void            CalcNormalization();
    double          UpdateLattice(double beta);
    void            PrintLattice(double beta);
    void            BetaRun(double beta);
    void            RunSimulation();
    void            InitSimulation();
    void            FreeSimulation();

    std::ofstream   outFile;
    Configuration   config;
    CRandomMersenne mersenne;

    int*            link;
    double*         scaleFactor;
    double*         sliceNormalization;
    double          normalization;

    double(*betaFunc)(double beta, int t, int tSize, double tScale, double(*scaleFactorFunc)(double, int tSize, double tScale));
    // betaFunc is a pointer to a function which returns a double, takes a (double, int, double, scaleFactoFunc), see below
    double(*scaleFactorFunc)(double, int tSize, double tScale);
    // scaleFactorFunc is a pointer to a function which returns a double, takes a (double, int, double)
};

#endif
