#ifndef GLOBALS_H
#define GLOBALS_H

typedef struct
{
    typedef enum {MINKOWSKI, RADIATION, MATTER, LAMBDA } cosmology_t;
    typedef enum {NORMAL, EFFECTIVE} mode_t;
    typedef enum {UNITY, RANDOM} init_t;

    int         t_size;
    int         s_size;
    cosmology_t cosmology;
    mode_t      mode;
    init_t      init;
    char        filename[1024];
    FILE*       fp;
    double      beta_min;
    double      beta_max;
    double      dbeta;
    int         iterations;
    double      t_scale;
} config_t;

typedef struct
{
    int*        link;
    double*     scale_factor;
    double*     slice_normalization;
    double      normalization;

    double(*beta_func)(double, int);
    double(*a_func)(double);
} simulation_t;

config_t        config;
simulation_t    simulation;

#define IDX(t, x, y, z, d) \
    (t * config.s_size * config.s_size * config.s_size * 4 + \
     x * config.s_size * config.s_size * 4 + \
     y * config.s_size * 4 + \
     z * 4 + \
     d)

CRandomMersenne mersenne(GetTickCount());

#endif
