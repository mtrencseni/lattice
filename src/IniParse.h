#ifndef INIPARSE_H
#define INIPARSE_H

static int handle_cosmology(config_t* config, const char* value)
{
#define MATCH(s) _stricmp(s, value) == 0

    if      (MATCH("minkowski"))    config->cosmology = config_t::MINKOWSKI;
    else if (MATCH("radiation"))    config->cosmology = config_t::RADIATION;
    else if (MATCH("matter"))       config->cosmology = config_t::MATTER;
    else if (MATCH("lambda"))       config->cosmology = config_t::LAMBDA;
    else                            return 0;

    return 1;

#undef MATCH
}

static int handle_mode(config_t* config, const char* value)
{
#define MATCH(s) _stricmp(s, value) == 0

    if      (MATCH("normal"))       config->mode = config_t::NORMAL;
    else if (MATCH("effective"))    config->mode = config_t::EFFECTIVE;
    else                            return 0;

    return 1;

#undef MATCH
}

static int handle_init(config_t* config, const char* value)
{
#define MATCH(s) _stricmp(s, value) == 0

    if      (MATCH("unity"))        config->init = config_t::UNITY;
    else if (MATCH("random"))       config->init = config_t::RANDOM;
    else                            return 0;

    return 1;

#undef MATCH
}

static int handler(void* user, const char* section, const char* name, const char* value)
{
#define MATCH(s) _stricmp(s, name) == 0

    config_t* config = (config_t*) user;

    if      (MATCH("t_size"))       config->t_size = atoi(value);
    else if (MATCH("s_size"))       config->s_size = atoi(value);
    else if (MATCH("cosmology"))    return handle_cosmology(config, value);
    else if (MATCH("mode"))         return handle_mode(config, value);
    else if (MATCH("init"))         return handle_init(config, value);
    else if (MATCH("output"))       strcpy(config->filename, value);
    else if (MATCH("beta_min"))     config->beta_min = atof(value);
    else if (MATCH("beta_max"))     config->beta_max = atof(value);
    else if (MATCH("dbeta"))        config->dbeta = atof(value);
    else if (MATCH("t_scale"))      config->t_scale = atof(value);
    else if (MATCH("iterations"))   config->iterations = atoi(value);
    else                            return 0;

    return 1;

#undef MATCH
}

#endif
