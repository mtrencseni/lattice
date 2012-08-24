#include <windows.h>
#include "simpleini/SimpleIni.h"
#include "Configuration.h"
#include "Lattice.h"

bool InitConfig(CSimpleIniA& ini, Configuration& config)
{
    config.tSize = ini.GetLongValue("", "t_size", -1);
    if (config.tSize == -1)
        return false;
    config.sSize = ini.GetLongValue("", "s_size", -1);
    if (config.sSize == -1)
        return false;
    
    std::string cosmology = ini.GetValue("", "cosmology");
    if      (cosmology == "minkowski")  config.cosmology = Configuration::MINKOWSKI;
    else if (cosmology == "radiation")  config.cosmology = Configuration::RADIATION;
    else if (cosmology == "matter")     config.cosmology = Configuration::MATTER;
    else if (cosmology == "lambda")     config.cosmology = Configuration::LAMBDA;
    else                                return false;

    std::string mode = ini.GetValue("", "mode");
    if      (mode == "normal")          config.mode = Configuration::NORMAL;
    else if (mode == "effective")       config.mode = Configuration::EFFECTIVE;
    else                                return false;

    std::string init = ini.GetValue("", "init");
    if      (init == "unity")                config.init = Configuration::UNITY;
    else if (init == "random")          config.init = Configuration::RANDOM;
    else                                return false;

    config.filename = ini.GetValue("", "output");
    if (config.filename == "")
        return false;
    config.betaMin = ini.GetDoubleValue("", "beta_min", -1);
    if (config.betaMin == -1)
        return false;
    config.betaMax = ini.GetDoubleValue("", "beta_max", -1);
    if (config.betaMax == -1)
        return false;
    config.dBeta = ini.GetDoubleValue("", "dbeta", -1);
    if (config.dBeta == -1)
        return false;
    config.tScale = ini.GetDoubleValue("", "t_scale", -1);
    if(config.tScale == -1)
        return false;
    config.iterations = ini.GetLongValue("", "iterations", -1);
    if (config.iterations == -1)
        return false;

    return true;
}

int main(int argc, char** argv)
{
    CSimpleIniA ini;
    Configuration config;
    Lattice lattice;

    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <config_file>" << std::endl;
        return 1;
    }

    ini.SetUnicode();
    if (ini.LoadFile(argv[1]) != SI_OK)
    {
        std::cout << "Can't load configuration from " << argv[1] << std::endl;
        return 1;
    }
    std::cout << "Configuration loaded from " << argv[1] << std::endl;

    if (!InitConfig(ini, config))
    {
        std::cout << "Can't load configuration from " << argv[1] << std::endl;
        return 1;
    }

    lattice.Run(config);

    return 0;
}
