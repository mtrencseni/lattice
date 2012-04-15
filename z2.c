/* Z_2 lattice gauge simulation			*/
/* Michael Creutz <creutz@bnl.gov>		*/
/* http://thy.phy.bnl.gov/~creutz/z2.c	*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// the lattice is of dimensions SIZE**4
#define T_SIZE      (20)
#define S_SIZE      (10)
#define ITER_MIN    (0)
#define ITER_MAX    (T_SIZE)
int link[T_SIZE][S_SIZE][S_SIZE][S_SIZE][4]; // last index gives link direction
double scale_factor[T_SIZE]; // a(t)
double normalization;

double a_radiation(double t) // a(t) ~ t^(1/2)
{
    if (t < T_SIZE/2)
        return sqrt(1.0 + t/16.0);
    else
        return sqrt(1.0 + (T_SIZE - t)/16.0);
}

// double a_matter(double t) // a(t) ~ t^(2/3)
// {
//     return pow(1 + t/8.0, 2.0/3.0);
// }
// 
// double a_lambda(double t) // a(t) ~ exp(t)
// {
//     return exp(1 + t/14);
// }

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

void init_unity()
{
	int d;
	int x[4];

	// set all links to unity

	for (x[0] = 0; x[0] < T_SIZE; x[0]++)
	{
		for (x[1] = 0; x[1] < S_SIZE; x[1]++)
		{
			for (x[2] = 0; x[2] < S_SIZE; x[2]++)
			{
				for (x[3] = 0; x[3] < S_SIZE; x[3]++)
				{
					for (d = 0; d < 4; d++)
						link[x[0]][x[1]][x[2]][x[3]][d] = 1;
				}
			}
		}
	}
}

void init_rand()
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
					    if (drand48() < 0.5)
						    link[x[0]][x[1]][x[2]][x[3]][d] = 1;
						else
						    link[x[0]][x[1]][x[2]][x[3]][d] = -1;
					}
				}
			}
		}
	}
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
	double action, staplesum;
    double a, gfactors;

    // do a Monte Carlo sweep; return energy

	normalization = 0.0;

	for (x[0] = ITER_MIN; x[0] < ITER_MAX; x[0]++)
	{
        a = scale_factor[x[0]]/scale_factor[0];
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
                            
                            normalization += 2 * gfactors;
                        }
                    }
                }
            }
        }
	}
}

double update(double beta)
{
	int	   x[4];
    int    u, v;
	double staple, staplesum;	
	double bplus, bminus, p;
	double action;
    double a, gfactors;

    // do a Monte Carlo sweep; return energy

	action = 0.0;

	for (x[0] = ITER_MIN; x[0] < ITER_MAX; x[0]++)
	{
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
                        bplus = exp(beta * staplesum);
                        bminus = 1 / bplus;
                        p = bplus / (bplus + bminus);
                        // printf("b+: %f\t", bplus);
                        // printf("b-: %f\t", bminus);
                        // printf("p : %f\n", p);
                        // the heatbath algorithm
                        if ( drand48() < p )
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
    
	//action /= (SIZE*SIZE*SIZE*SIZE*4*6);
    action /= normalization;
	return 1.0 - action;
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
    for (beta = 0.0; beta < max; beta += dbeta)
    {
        for (i = 0; i < num; i++)
            action = update(beta);
        printf("%g\t%f\n", beta, action); 
    }
}

void betarun(double beta, int num)
{
    double action;
    int i;

    for (i = 0; i < num; i++)
    {
        action = update(beta);
        printf("%i\t%f\n", i, action);
    }
}

void print_lattice()
{
    int x[4];
    int u, v;
	double bplus, bminus;
    int num_plus;
    int num_minus;
    double a, staple, staplesum, gfactors, action, sum_action;
    
    sum_action = 0.0;
    for (x[0] = ITER_MIN; x[0] < ITER_MAX; x[0]++)
	{
	    a = scale_factor[x[0]]/scale_factor[0];
        num_plus = 0;
        num_minus = 0;
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
                        if (link[x[0]][x[1]][x[2]][x[3]][u] == 1)
                            num_plus++;
                        else
                            num_minus++;
                            
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
        printf("t: %d\ta: %f\taction: %f\tsum_action: %f\t+: %d\t-: %d\n", x[0], a, action/normalization, sum_action/normalization, num_plus, num_minus);
    }
    printf("normalization: %f\n", normalization);
}

int main()
{
    srand48(1235L);

    //calc_scale_factor(a_minkowski);
    calc_scale_factor(a_radiation);
    calc_normalization();

    //init_unity();
    init_rand();
    heatcycle(1.0, 0.01, 1);

    //print_lattice();
    //betarun(0.2, 10);
    //print_lattice();
}
