#include "Lattice.h"
#ifdef WIN32
#include "windows.h"
#else
#include <time.h>
#endif

#define IDX(t, x, y, z, d) \
    (t * config.sSize * config.sSize * config.sSize * 4 + \
     x * config.sSize * config.sSize * 4 + \
     y * config.sSize * 4 + \
     z * 4 + \
     d)

Lattice::Lattice() : mersenne(
#ifdef WIN32
GetTickCount()
#else
clock_gettime(CLOCK_REALTIME)
#endif
) { }

void Lattice::MoveUp(int x[], int d)
{
    int size;
    
    if (d == 0)
        size = config.tSize;
    else
        size = config.sSize;
    
    x[d] += 1;
    
    if (x[d] >= size)
        x[d] -= size; 
}

void Lattice::MoveDown(int x[], int d)
{
    int size;
    
    if (d == 0)
        size = config.tSize;
    else
        size = config.sSize;

    x[d] -= 1;
    
    if (x[d] < 0)
        x[d] += size;
}

void Lattice::InitRand(double r)
{
    int x[4], d, v;

    // randomize all links

    for (x[0] = 0; x[0] < config.tSize; x[0]++)
    {
        for (x[1] = 0; x[1] < config.sSize; x[1]++)
        {
            for (x[2] = 0; x[2] < config.sSize; x[2]++)
            {
                for (x[3] = 0; x[3] < config.sSize; x[3]++)
                {
                    for (d = 0; d < 4; d++)
                    {
                        if (mersenne.Random() < r)
                            v = 1;
                        else
                            v = -1;

                        link[IDX(x[0], x[1], x[2], x[3], d)] = v;
                    }
                }
            }
        }
    }
}

void Lattice::InitLattice()
{
    if (config.init == Configuration::UNITY)
        InitRand(1.0); // sets all links to 1
    else
        InitRand(0.5); // sets links to 0 or 1 with 1/2 probability
}

void Lattice::CalcScaleFactor()
{
    int     t;
    double(*a)(double, int tSize, double tScale);

    if (config.mode == Configuration::EFFECTIVE)
        a = aMinkowski;
    else
        a = scaleFactorFunc;
    
    for (t = 0; t < config.tSize; t++)
        scaleFactor[t] = a(t, config.tSize, config.tScale);
}

void Lattice::CalcNormalization()
{
    int	   x[4], u, v;
    double sliceAction, a, gFactors;

   normalization = 0.0;

    for (x[0] = 0; x[0] < config.tSize; x[0]++)
    {
        a = scaleFactor[x[0]]/scaleFactor[0];
        sliceAction = 0.0;

        for (x[1] = 0; x[1] < config.sSize; x[1]++)
        {
            for (x[2] = 0; x[2] < config.sSize; x[2]++)
            {
                for (x[3] = 0; x[3] < config.sSize; x[3]++)
                {
                    for (u = 0; u < 4; u++)
                    {
                        for (v = 0; v < 4; v++)
                        {
                            if (u == v)
                                continue;
                                
                            if (u == 0 || v == 0)
                                gFactors = a*a*a*a*a;       // a(t)^5
                            else
                                gFactors = a*a*a*a*a*a*a;   // a(t)^7
                            
                            sliceAction += 2 * gFactors;
                        }
                    }
                }
            }
        }
        
        sliceNormalization[x[0]] = sliceAction;
        normalization += sliceAction;
    }
}

double Lattice::UpdateLattice(double beta)
{
    int	   x[4], u, v;
    double a, p, action, bPlus, bMinus, betaSlice, gFactors, staple, stapleSum;

    // do a Monte Carlo sweep; return energy

    action = 0.0;

    for (x[0] = 0; x[0] < config.tSize; x[0]++)
    {
        betaSlice = betaFunc(beta, x[0], config.tSize, config.tScale, scaleFactorFunc);

        a = scaleFactor[x[0]]/scaleFactor[0];
        for (x[1] = 0; x[1] < config.sSize; x[1]++)
        {
            for (x[2] = 0; x[2] < config.sSize; x[2]++)
            {
                for (x[3] = 0; x[3] < config.sSize; x[3]++)
                {
                    for (u = 0; u < 4; u++)
                    {
                        stapleSum = 0.0;
                        for (v = 0; v < 4; v++)
                        {
                            if (u == v)
                                continue;
                                
                            if (u == 0 || v == 0)
                                gFactors = a*a*a*a*a;       // a(t)^5
                            else
                                gFactors = a*a*a*a*a*a*a;   // a(t)^7

                            /*  move around:
                                            6--5
                                ^            |  |
                                | v          1--4
                                |            |  |
                                -----> u     2--3  */
                                 
                            // plaquette 1234
                            staple = 1;
                            MoveDown(x, v);
                            staple *= link[IDX(x[0], x[1], x[2], x[3], u)]; /* 23 */
                            staple *= link[IDX(x[0], x[1], x[2], x[3], v)]; /* 12 */
                            MoveUp(x, u);
                            staple *= link[IDX(x[0], x[1], x[2], x[3], v)]; /* 34 */
                            stapleSum += staple * gFactors;

                            MoveUp(x, v);
                                                        
                            // plaquette 1456
                            staple = 1;
                            staple *= link[IDX(x[0], x[1], x[2], x[3], v)];  /* 45 */
                            MoveUp(x, v);
                            MoveDown(x, u);
                            staple *= link[IDX(x[0], x[1], x[2], x[3], u)]; /* 56 */
                            MoveDown(x, v);
                            staple *= link[IDX(x[0], x[1], x[2], x[3], v)]; /* 61 */

                            stapleSum += staple * gFactors;
                        }
                    
                        // calculate the Boltzmann weight
                        bPlus = exp(betaSlice * stapleSum);
                        bMinus = 1 / bPlus;
                        p = bPlus / (bPlus + bMinus);
                        // the heatbath algorithm
                        if (mersenne.Random() < p)
                        {
                            link[IDX(x[0], x[1], x[2], x[3], u)] = 1;
                            action += stapleSum;
                        }
                        else
                        {
                            link[IDX(x[0], x[1], x[2], x[3], u)] = -1;
                            action -= stapleSum;
                        }
                    }
                }
            }
        }
    }
    
    action /= normalization;
    return 1.0 - action;
}


void Lattice::PrintLattice(double beta)
{
    int     x[4], u, v;
    double  a, staple, stapleSum, gFactors, action;

    for (x[0] = 0; x[0] < config.tSize; x[0]++)
    {
        a = scaleFactor[x[0]]/scaleFactor[0];
        action = 0.0;
        for (x[1] = 0; x[1] < config.sSize; x[1]++)
        {
            for (x[2] = 0; x[2] < config.sSize; x[2]++)
            {
                for (x[3] = 0; x[3] < config.sSize; x[3]++)
                {
                    for (u = 0; u < 4; u++)
                    {
                        stapleSum = 0.0;
                        for (v = 0; v < 4; v++)
                        {
                            if (u == v)
                                continue;
                                
                            if (u == 0 || v == 0)
                                gFactors = a*a*a*a*a;       // a(t)^5
                            else
                                gFactors = a*a*a*a*a*a*a;   // a(t)^7

                            /*  move around:
                                           6--5
                              ^            |  |
                              | v          1--4
                              |            |  |
                              -----> u     2--3  */
                                 
                            // plaquette 1234
                            staple = 1; //link[x[0]][x[1]][x[2]][x[3]][u];  /* 41 */
                            MoveDown(x, v);
                            staple *= link[IDX(x[0], x[1], x[2], x[3], u)]; /* 23 */
                            staple *= link[IDX(x[0], x[1], x[2], x[3], v)]; /* 12 */
                            MoveUp(x, u);
                            staple *= link[IDX(x[0], x[1], x[2], x[3], v)]; /* 34 */
                            stapleSum += staple * gFactors;
                            
                            MoveUp(x, v);
                                                        
                            // plaquette 1456
                            staple = 1;
                            staple *= link[IDX(x[0], x[1], x[2], x[3], v)];  /* 45 */
                            MoveUp(x, v);
                            MoveDown(x, u);
                            staple *= link[IDX(x[0], x[1], x[2], x[3], u)]; /* 56 */
                            MoveDown(x, v);
                            staple *= link[IDX(x[0], x[1], x[2], x[3], v)]; /* 61 */
                            //staple *= link[x[0]][x[1]][x[2]][x[3]][u]; /* 14 */
                            stapleSum += staple * gFactors;
                        }
                        stapleSum *= link[IDX(x[0], x[1], x[2], x[3], u)];
                        action += stapleSum;
                    }
                }
            }
        }
        
        //printf(" %d", x[0]);
        std::cout << " " << x[0];

        if (x[0] < config.tSize / 2)
            outFile << x[0] << " "  << beta << " " << 1.0 - action/sliceNormalization[x[0]] << std::endl;
            //fprintf(config.fp, "%d %f %f\n", x[0], beta, 1.0 - action/sliceNormalization[x[0]]);
    }
}

void Lattice::BetaRun(double beta)
{
    int     i;
    double  action;

    for (i = 0; i < config.iterations; i++)
        action = UpdateLattice(beta);
}

void Lattice::RunSimulation()
{
    double beta;
    
    for (beta = config.betaMin; beta <= config.betaMax; beta += config.dBeta)
    {
        //printf("\nDoing %s, beta = %0.3f, t = ", config.filename, beta);
        std::cout << std::endl << "Doing " << config.filename << ", beta = " << beta << ", t = ";
        InitLattice();
        BetaRun(beta);
        PrintLattice(beta);
    }
}

void Lattice::InitSimulation()
{
    link = (int*) malloc(sizeof(int) * 4 * config.tSize * config.sSize * config.sSize * config.sSize);
    scaleFactor = (double*) malloc(sizeof(double) * config.tSize);
    sliceNormalization = (double*) malloc(sizeof(double) * config.tSize);
    
    if (config.mode == Configuration::NORMAL)
        betaFunc = betaNormal;
    else
        betaFunc = betaEff;

    if (config.cosmology == Configuration::MINKOWSKI)
        scaleFactorFunc = aMinkowski;
    else if (config.cosmology == Configuration::RADIATION)
        scaleFactorFunc = aRadiation;
    else if (config.cosmology == Configuration::MATTER)
        scaleFactorFunc = aMatter;
    else
        scaleFactorFunc = aLambda;

    CalcScaleFactor();
    CalcNormalization();
}

void Lattice::FreeSimulation()
{
    free(link);
    free(scaleFactor);
    free(sliceNormalization);
}

void Lattice::Run(Configuration& config_)
{
    config = config_;
 
    outFile.open(config.filename, std::ios::out);
    if(!outFile.is_open())
    {
        std::cout << "Unable to open file " << config.filename << " for writing" << std::endl;
        return;
    }
    
    InitSimulation();
    RunSimulation();
    FreeSimulation();

    outFile.close();
}
