#include <windows.h>
#include "mersenne/randomc.h"
#include "ini.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Globals.h"
#include "IniParse.h"

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
    a = simulation.a_func(t)/simulation.a_func(0);
    return beta * (a*a*a*a*a + a*a*a*a*a*a*a) / 2.0;
}

double a_radiation(double t) // a(t) ~ t^(1/2)
{
    if (t < config.t_size/2)
       return sqrt(1.0 + t/config.t_scale);
    else
        return sqrt(1.0 + (config.t_size - t)/config.t_scale);
}

double a_matter(double t) // a(t) ~ t^(2/3)
{
    if (t < config.t_size/2)
        return pow(1 + t/config.t_scale, 2.0/3.0);
    else
        return pow(1 + (config.t_size - t)/config.t_scale, 2.0/3.0);
}

double a_lambda(double t) // a(t) ~ exp(t)
{
    if (t < config.t_size/2)
        return exp(1 + t/config.t_scale);
    else
        return exp(1 + (config.t_size - t)/config.t_scale);
        
}

double a_minkowski(double t)
{
    return 1.0;
}

void moveup(int x[], int d)
{
    int size;
    
    if (d == 0)
        size = config.t_size;
    else
        size = config.s_size;
    
    x[d] += 1;
    
    if (x[d] >= size)
        x[d] -= size; 
}

void movedown(int x[], int d)
{
    int size;
    
    if (d == 0)
        size = config.t_size;
    else
        size = config.s_size;

    x[d] -= 1;
    
    if (x[d] < 0)
        x[d] += size;
}

void init_rand(double r)
{
    int d;
    int x[4];
    int v;

    // randomize all links

    for (x[0] = 0; x[0] < config.t_size; x[0]++)
    {
        for (x[1] = 0; x[1] < config.s_size; x[1]++)
        {
            for (x[2] = 0; x[2] < config.s_size; x[2]++)
            {
                for (x[3] = 0; x[3] < config.s_size; x[3]++)
                {
                    for (d = 0; d < 4; d++)
                    {
                        if (rnd() < r)
                            v = 1;
                        else
                            v = -1;

                        simulation.link[IDX(x[0], x[1], x[2], x[3], d)] = v;
                    }
                }
            }
        }
    }
}

void init_lattice()
{
    if (config.init == config_t::UNITY)
        init_rand(1.0); // sets all links to 1
    else
        init_rand(0.5); // sets links to 0 or 1 with 1/2 probability
}

void calc_scale_factor()
{
    int     t;
    double  (*a_func)(double);

    if (config.mode == config_t::EFFECTIVE)
        a_func = a_minkowski;
    else
        a_func = simulation.a_func;
    
    for (t = 0; t < config.t_size; t++)
        simulation.scale_factor[t] = a_func(t);
}

void calc_normalization()
{
    int	   x[4];
    int    u, v;
    double slice_action;
    double a, gfactors;

   simulation.normalization = 0.0;

    for (x[0] = 0; x[0] < config.t_size; x[0]++)
    {
        a = simulation.scale_factor[x[0]]/simulation.scale_factor[0];
        slice_action = 0.0;

        for (x[1] = 0; x[1] < config.s_size; x[1]++)
        {
            for (x[2] = 0; x[2] < config.s_size; x[2]++)
            {
                for (x[3] = 0; x[3] < config.s_size; x[3]++)
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
        
        simulation.slice_normalization[x[0]] = slice_action;
        simulation.normalization += slice_action;
    }
}

double update(double beta)
{
    int	   x[4];
    int    u, v;
    double beta_slice;
    double staple, staplesum;	
    double bplus, bminus, p;
    double action;
    double a, gfactors;

    // do a Monte Carlo sweep; return energy

    action = 0.0;

    for (x[0] = 0; x[0] < config.t_size; x[0]++)
    {
        beta_slice = simulation.beta_func(beta, x[0]);

        a = simulation.scale_factor[x[0]]/simulation.scale_factor[0];
        for (x[1] = 0; x[1] < config.s_size; x[1]++)
        {
            for (x[2] = 0; x[2] < config.s_size; x[2]++)
            {
                for (x[3] = 0; x[3] < config.s_size; x[3]++)
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
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], u)]; /* 23 */
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], v)]; /* 12 */
                            moveup(x, u);
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], v)]; /* 34 */
                            staplesum += staple * gfactors;

                            moveup(x, v);
                                                        
                            // plaquette 1456
                            staple = simulation.link[IDX(x[0], x[1], x[2], x[3], v)];  /* 45 */
                            moveup(x, v);
                            movedown(x, u);
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], u)]; /* 56 */
                            movedown(x, v);
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], v)]; /* 61 */

                            staplesum += staple * gfactors;
                        }
                    
                        // calculate the Boltzmann weight
                        bplus = exp(beta_slice * staplesum);
                        bminus = 1 / bplus;
                        p = bplus / (bplus + bminus);
                        // the heatbath algorithm
                        if ( rnd() < p )
                        {
                            simulation.link[IDX(x[0], x[1], x[2], x[3], u)] = 1;
                            action += staplesum;
                        }
                        else
                        {
                            simulation.link[IDX(x[0], x[1], x[2], x[3], u)] = -1;
                            action -= staplesum;
                        }
                    }
                }
            }
        }
    }
    
    action /= simulation.normalization;
    return 1.0 - action;
}

void print_lattice(double beta)
{
    int     x[4];
    int     u, v;
    double  a, staple, staplesum, gfactors, action;

    for (x[0] = 0; x[0] < config.t_size; x[0]++)
    {
        a = simulation.scale_factor[x[0]]/simulation.scale_factor[0];
        action = 0.0;
        for (x[1] = 0; x[1] < config.s_size; x[1]++)
        {
            for (x[2] = 0; x[2] < config.s_size; x[2]++)
            {
                for (x[3] = 0; x[3] < config.s_size; x[3]++)
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
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], u)]; /* 23 */
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], v)]; /* 12 */
                            moveup(x, u);
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], v)]; /* 34 */
                            staplesum += staple * gfactors;
                            
                            moveup(x, v);
                                                        
                            // plaquette 1456
                            staple = simulation.link[IDX(x[0], x[1], x[2], x[3], v)];  /* 45 */
                            moveup(x, v);
                            movedown(x, u);
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], u)]; /* 56 */
                            movedown(x, v);
                            staple *= simulation.link[IDX(x[0], x[1], x[2], x[3], v)]; /* 61 */
                            //staple *= link[x[0]][x[1]][x[2]][x[3]][u]; /* 14 */
                            staplesum += staple * gfactors;
                        }
                        staplesum *= simulation.link[IDX(x[0], x[1], x[2], x[3], u)];
                        action += staplesum;
                    }
                }
            }
        }
        
        printf(" %d", x[0]);

        if (x[0] < config.t_size / 2)
            fprintf(config.fp, "%d %f %f\n", x[0], beta, 1.0 - action/simulation.slice_normalization[x[0]]);
    }
}

void betarun(double beta)
{
    int     i;
    double  action;

    for (i = 0; i < config.iterations; i++)
        action = update(beta);
}

void run_simulation()
{
    double beta;
    
    for (beta = config.beta_min; beta <= config.beta_max; beta += config.dbeta)
    {
        printf("\nDoing %s, beta = %0.3f, t = ", config.filename, beta);
        init_lattice();
        betarun(beta);
        print_lattice(beta);
    }
}

void init_simulation()
{
    simulation.link = (int*) malloc(sizeof(int) * 4 * config.t_size * config.s_size * config.s_size * config.s_size);
    simulation.scale_factor = (double*) malloc(sizeof(double) * config.t_size);
    simulation.slice_normalization = (double*) malloc(sizeof(double) * config.t_size);
    
    if (config.mode == config_t::NORMAL)
        simulation.beta_func = beta_normal;
    else
        simulation.beta_func = beta_eff;

    if (config.cosmology == config_t::MINKOWSKI)
        simulation.a_func = a_minkowski;
    else if (config.cosmology == config_t::RADIATION)
        simulation.a_func = a_radiation;
    else if (config.cosmology == config_t::MATTER)
        simulation.a_func = a_matter;
    else
        simulation.a_func = a_lambda;

    calc_scale_factor();
    calc_normalization();
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        printf("Usage: %s <config_file>", argv[0]);
        return 1;
    }

    if (ini_parse(argv[1], handler, &config) < 0)
    {
        printf("Can't load configuration from argv[1]\n");
        return 1;
    }
    
    printf("Configuration loaded from %s\n", argv[1]);

    config.fp = fopen(config.filename, "w+");
    if (config.fp == NULL)
    {
        printf("Unable to open file %s for writing\n", config.filename);
        return 1;
    }

    init_simulation();
    run_simulation();

    fclose(config.fp);

    return 0;
}

/*
gnuplot heatmap from x y z lines:

set tics out
set xrange[0:20]
set yrange[0:0.5]
set xlabel "t"
set ylabel "$\\beta$"
set cbrange[0:1]
set palette defined (0 "white", 1 "#303030")

set term windows
set title "(a) radiation"
plot "d:/lattice/data/rad.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/(sqrt(1+x/32)**5+sqrt(1+x/32)**7) linecolor rgb "black" lw 6 title ""
set term epslatex
set output "d:/lattice/data/rad.tex"
replot
set output

set term windows
set title "(b) matter"
plot "d:/lattice/data/mat.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/(((1+x/32)**(0.66))**5+((1+x/32)**(0.66))**7) linecolor rgb "black" lw 6 title ""
set term epslatex
set output "d:/lattice/data/mat.tex"
replot
set output

set term windows
set title "(c) lambda"
plot "d:/lattice/data/lam.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/((exp(1+x/32)/exp(1))**5+(exp(1+x/32)/exp(1))**7) lc rgb "black" lw 6 title ""
set term epslatex
set output "d:/lattice/data/lam.tex"
replot
set output

DIFFERENCE PLOTS:

set palette defined (0 "white", 0.05 "#303030")
set cbrange[0:0.05]

set term windows
set title "(a) radiation"
plot "d:/lattice/data/diffrad.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/(sqrt(1+x/32)**5+sqrt(1+x/32)**7) linecolor rgb "black" lw 6 title ""
set term epslatex color
set output "d:/lattice/data/diffrad.tex"
replot
set output

set term windows
set title "(b) matter"
plot "d:/lattice/data/diffmat.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/(((1+x/32)**(0.66))**5+((1+x/32)**(0.66))**7) linecolor rgb "black" lw 6 title ""
set term epslatex color
set output "d:/lattice/data/diffmat.tex"
replot
set output

set term windows
set title "(c) lambda"
plot "d:/lattice/data/difflam.dat" u ($1+0.5):($2+0.0125):($3) w image title "", 0.88/((exp(1+x/32)/exp(1))**5+(exp(1+x/32)/exp(1))**7) linecolor rgb "black" lw 6 title ""
set term epslatex color
set output "d:/lattice/data/difflam.tex"
replot
set output

FIXED BETA:

set tics out
set xlabel "t"
set ylabel "$E$"
set term windows
set title ""
set xrange[0:10]
set yrange[0:0.8]
unset arrow
set arrow from 6.57,0.8 to 6.57,0 nohead lc rgb "black"
plot "d:/lattice/data/rad_beta025.dat" u ($1)/10:($3) linecolor rgb "black" pt 1 title "", "d:/lattice/data/radeff_beta025.dat" u ($1)/10:($3) linecolor rgb "black" pt 4 title ""
set term epslatex
set output "d:/lattice/data/rad_beta025.tex"
replot
set output


set tics out
set xlabel "t"
set ylabel "$E$"
set term windows
set title ""
set xrange[5:15]
set yrange[0:0.8]
unset arrow
set arrow from 6.57,0.8 to 6.57,0 nohead lc rgb "black"
set arrow from (40-6.57),0.8 to (40-6.57),0.3 nohead lc rgb "black"
plot "d:/lattice/data/rad_beta025.dat" u ($1)/10:($3) linecolor rgb "black" pt 1 title "", "d:/lattice/data/rad_beta025.dat" u (400-$1)/10:($3) linecolor rgb "gray" pt 6 title ""
set term epslatex
set output "d:/lattice/data/rad_beta025_hyst.tex"
replot
set output




set tics out
set xlabel "$\beta$"
set ylabel "$E$"
set term windows
set title ""
set xrange[0:2]
set yrange[0:1]
unset arrow
set arrow from 0.4,1 to 0.4,0 nohead lc rgb "black"
plot "d:/lattice/data/tbeta2.dat" u ($1)/40:($3) linecolor rgb "black" pt 1 title "", "d:/lattice/data/tbeta2.dat" u (80-$1)/40:($3) linecolor rgb "gray" pt 6 title ""



set arrow from 0.44,0 to 0.44,1 nohead lc rgb "black"

*/
