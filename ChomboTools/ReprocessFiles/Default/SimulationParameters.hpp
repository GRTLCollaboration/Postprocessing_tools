/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    // For the Interpolator test we don't need many parameters
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
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

        // extraction params
        dx.fill(coarsest_dx);
        origin.fill(coarsest_dx / 2.0);
    }

    int num_files, start_file, end_file;
    std::array<double, CH_SPACEDIM> origin,
        dx; // location of coarsest origin and dx
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
