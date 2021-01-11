/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif

// General includes:
#include "parstream.H" //Gives us pout()
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/time.h>
using std::endl;

// Problem specific includes:
#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "ReprocessingLevel.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
#include "UserVariables.hpp"

int runReprocessingTool(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);
    SimulationParameters sim_params(pp);

    // Setup the initial object (from restart_file checkpoint)
    BHAMR bh_amr;
    DefaultLevelFactory<ReprocessingLevel> empty_level_fact(bh_amr, sim_params);
    setupAMRObject(bh_amr, empty_level_fact);

    // set up interpolator
    AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    bh_amr.set_interpolator(&interpolator);

#if USE_AHFINDER
    AHSphericalGeometry sph(sim_params.center);
    bh_amr.m_ah_finder.add_ah(sph, sim_params.AH_initial_guess,
                              sim_params.AH_params, "stats", "coords", false);
#endif

    // now loop over chk files
    for (int ifile = sim_params.start_file; ifile <= sim_params.end_file;
         ifile += sim_params.checkpoint_interval)
    {
        // set up the file from next checkpoint
        std::ostringstream current_file;
        current_file << std::setw(6) << std::setfill('0') << ifile;
        std::string restart_file(sim_params.checkpoint_prefix +
                                 current_file.str() + ".3d.hdf5");
        // std::cout << "Opening file '" << restart_file << "'" << std::endl;
        HDF5Handle handle(restart_file, HDF5Handle::OPEN_RDONLY);
        bh_amr.setupForRestart(handle);
        handle.close();
    }

    bh_amr.conclude();

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runReprocessingTool(argc, argv);

    mainFinalize();
    return status;
}
