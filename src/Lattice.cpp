/* Z_2 lattice gauge simulation			*/
/* Michael Creutz <creutz@bnl.gov>		*/
/* http://thy.phy.bnl.gov/~creutz/z2.c	*/

#include <windows.h>
#include "mersenne/randomc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define T_SIZE      (30)
#define S_SIZE      (20)
#define ITER_MIN    (0)
#define ITER_MAX    (T_SIZE)
int link[T_SIZE][S_SIZE][S_SIZE][S_SIZE][4]; // last index gives link direction
double scale_factor[T_SIZE]; // a(t)
double slice_normalization[T_SIZE];
double normalization;
double(*beta_func)(double, int);
double(*a_func_eff)(double);

CRandomMersenne mersenne(GetTickCount());

double rnd()
{
    return mersenne.Random();
}

double beta_normal(double beta, int t)
{
    return beta;
}

double beta_eff(double beta, int t)
{
    double a;
    a = a_func_eff(t)/a_func_eff(0);
    return beta * (a*a*a*a*a + a*a*a*a*a*a*a) / 2.0;
    //return beta * (a*a*a*a*a*a);
}

double a_radiation(double t) // a(t) ~ t^(1/2)
{
    // basic plot
    /*if (t < T_SIZE/2)
       return sqrt(1.0 + t/16.0);
    else
        return sqrt(1.0 + (T_SIZE - t)/16.0);*/

    // high res plot
    return sqrt(1.0 + t/(8*16.0));
}

double a_matter(double t) // a(t) ~ t^(2/3)
{
    if (t < T_SIZE/2)
        return pow(1 + t/16.0, 2.0/3.0);
    else
        return pow(1 + (T_SIZE - t)/16.0, 2.0/3.0);
}

double a_lambda(double t) // a(t) ~ exp(t)
{
    if (t < T_SIZE/2)
        return exp(1 + t/16.0);
    else
        return exp(1 + (T_SIZE - t)/16.0);
        
}

double a_minkowski(double t)
{
    return 1.0;
}

void moveup(int x[], int d)
{
    int size;
    
    if (d == 0)
        size = T_SIZE;
    else
        size = S_SIZE;
    
    x[d] += 1;
    
    if (x[d] >= size)
        x[d] -= size; 
}

void movedown(int x[],int d)
{
    int size;
    
    if (d == 0)
        size = T_SIZE;
    else
        size = S_SIZE;

    x[d] -= 1;
    
    if (x[d] < 0)
        x[d] += size;
}

void init_rand(double r)
{
    int d;
    int x[4];

    // randomize all links

    for (x[0] = 0; x[0] < T_SIZE; x[0]++)
    {
        for (x[1] = 0; x[1] < S_SIZE; x[1]++)
        {
            for (x[2] = 0; x[2] < S_SIZE; x[2]++)
            {
                for (x[3] = 0; x[3] < S_SIZE; x[3]++)
                {
                    for (d = 0; d < 4; d++)
                    {
                        if (rnd() < r)
                            link[x[0]][x[1]][x[2]][x[3]][d] = 1;
                        else
                            link[x[0]][x[1]][x[2]][x[3]][d] = -1;
                    }
                }
            }
        }
    }
}

void init_unity()
{
    init_rand(1.0); // sets all links to 1
}

void init_uniform()
{
    init_rand(0.5); // sets links to 0 or 1 with 1/2 probability
}

void calc_scale_factor( double(*a_func)(double) )
{
    int t;
    
    for (t = 0; t < T_SIZE; t++)
        scale_factor[t] = a_func(t);
}

void calc_normalization()
{
    int	   x[4];
    int    u, v;
    double slice_action;
    double a, gfactors;

   normalization = 0.0;

    for (x[0] = ITER_MIN; x[0] < ITER_MAX; x[0]++)
    {
        a = scale_factor[x[0]]/scale_factor[0];
        slice_action = 0.0;
        for (x[1] = 0; x[1] < S_SIZE; x[1]++)
        {
            for (x[2] = 0; x[2] < S_SIZE; x[2]++)
            {
                for (x[3] = 0; x[3] < S_SIZE; x[3]++)
                {
                    for (u = 0; u < 4; u++)
                    {
                        for (v = 0; v < 4; v++)
                        {
                            if (u == v)
                                continue;
                                
                            if (u == 0 || v == 0)
                                gfactors = a*a*a*a*a;       // a(t)^5
                            else
                                gfactors = a*a*a*a*a*a*a;   // a(t)^7
                            
                            slice_action += 2 * gfactors;
                        }
                    }
                }
            }
        }
        slice_normalization[x[0]] = slice_action;
        normalization += slice_action;
    }
}

double update(double beta)
{
    int	   x[4];
    int    u, v;
    double beta_eff;
    double staple, staplesum;	
    double bplus, bminus, p;
    double action;
    double a, gfactors;

    // do a Monte Carlo sweep; return energy

    beta_eff = beta;
    action = 0.0;

    for (x[0] = ITER_MIN; x[0] < ITER_MAX; x[0]++)
    {
        a = a_matter(x[0]) / a_matter(0);
        beta_eff = beta * (a*a*a*a*a + a*a*a*a*a*a*a) / 2.0;
        a = scale_factor[x[0]]/scale_factor[0];
        for (x[1] = 0; x[1] < S_SIZE; x[1]++)
        {
            for (x[2] = 0; x[2] < S_SIZE; x[2]++)
            {
                for (x[3] = 0; x[3] < S_SIZE; x[3]++)
                {
                    for (u = 0; u < 4; u++)
                    {
                        staplesum = 0.0;
                        for (v = 0; v < 4; v++)
                        {
                            if (u == v)
                                continue;
                                
                            if (u == 0 || v == 0)
                                gfactors = a*a*a*a*a;       // a(t)^5
                            else
                                gfactors = a*a*a*a*a*a*a;   // a(t)^7

                            /*  move around:
                                           6--5
                              ^            |  |
                              | v          1--4
                              |            |  |
                              -----> u     2--3  */
                                 
                            // plaquette 1234
                            staple = 1;
                            movedown(x, v);
                            staple *= link[x[0]][x[1]][x[2]][x[3]][u]; /* 23 */
                            staple *= link[x[0]][x[1]][x[2]][x[3]][v]; /* 12 */
                            moveup(x, u);
                            staple *= link[x[0]][x[1]][x[2]][x[3]][v]; /* 34 */
                            staplesum += staple * gfactors;

                            moveup(x, v);
                                                        
                            // plaquette 1456
                            staple = link[x[0]][x[1]][x[2]][x[3]][v];  /* 45 */
                            moveup(x, v);
                            movedown(x, u);
                            staple *= link[x[0]][x[1]][x[2]][x[3]][u]; /* 56 */
                            movedown(x, v);
                            staple *= link[x[0]][x[1]][x[2]][x[3]][v]; /* 61 */

                            staplesum += staple * gfactors;
                        }
                    
                        // calculate the Boltzmann weight
                        bplus = exp(beta_func(beta, x[0]) * staplesum);
                        bminus = 1 / bplus;
                        p = bplus / (bplus + bminus);
                        // the heatbath algorithm
                        if ( rnd() < p )
                        {
                            link[x[0]][x[1]][x[2]][x[3]][u] = 1;
                            action += staplesum;
                        }
                        else
                        {
                            link[x[0]][x[1]][x[2]][x[3]][u] = -1;
                            action -= staplesum;
                        }
                    }
                }
            }
        }
    }
    
    action /= normalization;
    return action;
}

void print_lattice(double beta)
{
    int x[4];
    int u, v;
    double a, staple, staplesum, gfactors, action, sum_action;
    sum_action = 0.0;
    for (x[0] = ITER_MIN; x[0] < ITER_MAX; x[0]++)
    {
        a = scale_factor[x[0]]/scale_factor[0];
        action = 0.0;
        for (x[1] = 0; x[1] < S_SIZE; x[1]++)
        {
            for (x[2] = 0; x[2] < S_SIZE; x[2]++)
            {
                for (x[3] = 0; x[3] < S_SIZE; x[3]++)
                {
                    for (u = 0; u < 4; u++)
                    {
                        staplesum = 0.0;
                            
                        for (v = 0; v < 4; v++)
                        {
                            if (u == v)
                                continue;
                                
                            if (u == 0 || v == 0)
                                gfactors = a*a*a*a*a;       // a(t)^5
                            else
                                gfactors = a*a*a*a*a*a*a;   // a(t)^7

                            /*  move around:
                                           6--5
                              ^            |  |
                              | v          1--4
                              |            |  |
                              -----> u     2--3  */
                                 
                            // plaquette 1234
                            staple = 1; //link[x[0]][x[1]][x[2]][x[3]][u];  /* 41 */
                            movedown(x, v);
                            staple *= link[x[0]][x[1]][x[2]][x[3]][u]; /* 23 */
                            staple *= link[x[0]][x[1]][x[2]][x[3]][v]; /* 12 */
                            moveup(x, u);
                            staple *= link[x[0]][x[1]][x[2]][x[3]][v]; /* 34 */
                            staplesum += staple * gfactors;
                            
                            moveup(x, v);
                                                        
                            // plaquette 1456
                            staple = link[x[0]][x[1]][x[2]][x[3]][v];  /* 45 */
                            moveup(x, v);
                            movedown(x, u);
                            staple *= link[x[0]][x[1]][x[2]][x[3]][u]; /* 56 */
                            movedown(x, v);
                            staple *= link[x[0]][x[1]][x[2]][x[3]][v]; /* 61 */
                            //staple *= link[x[0]][x[1]][x[2]][x[3]][u]; /* 14 */
                            staplesum += staple * gfactors;
                        }
                        action += fabs(staplesum);
                    }
                }
            }
        }
        sum_action += action;
        printf("%d %f %f\n",
         x[0], beta, action/slice_normalization[x[0]]);
    }
}

void heatcycle(double max, double dbeta, int num)
{
    double beta, action;
    int i;

    // heat it up
    for (beta = max; beta > 0.0; beta -= dbeta)
    {
        for (i = 0; i < num; i++)
            action = update(beta);
        printf("%g\t%f\n", beta, action); 
    }
    
    // cool it down
    // for (beta = 0.0; beta < max; beta += dbeta)
    // {
    //     for (i = 0; i < num; i++)
    //         action = update(beta);
    //     printf("%g\t%f\n", beta, action); 
    // }
}

void betarun(double beta, int num)
{
    int     i;
    double  action;

    for (i = 0; i < num; i++)
    {
        action = update(beta);
        //printf("%i\t%f\n", i, action);
    }    
}

void betaruns(double min, double max, double dbeta, int num, void(*lattice_init_func)())
{
    double beta;
    
    for (beta = min; beta <= max; beta += dbeta)
    {
        lattice_init_func();
        betarun(beta, num);
        print_lattice(beta);
    }
}

int main()
{
    calc_scale_factor(a_minkowski);
    calc_normalization();
    a_func_eff = a_radiation;
    beta_func = beta_eff;

    //heatcycle(1.0, 0.01, 10);

    //betaruns(0.15, 0.15, 1, 10, init_unity);
    //betaruns(0, 0.5, 0.025, 10, init_unity);
    betaruns(0.25, 0.35, 0.0025, 10, init_unity); // high res plot for radiation case
}

/*
gnuplot heatmap from x y z lines:

set term x11
set title "Normalized action per spacelike slice for different betas"
set xlabel "Time"
set ylabel "Beta"
set palette defined (0 "blue", 1 "red")
splot "3d.dat"
plot "3dmat.dat" u ($1+0.5):($2+0.0125):($3) w image
plot "3drad_hires.dat" u ($1+0.5):($2+0.00125):($3) w image

*/