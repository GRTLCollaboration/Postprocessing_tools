/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    // For the Interpolator test we don't need many parameters
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // Grid setup
        pp.load("num_files", num_files, -1);
        pp.load("start_file", start_file);

        if (num_files == -1)
            pp.load("end_file", end_file);
        else
            end_file = start_file + num_files * checkpoint_interval;

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess);
#endif
    }

    int num_files, start_file, end_file;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
