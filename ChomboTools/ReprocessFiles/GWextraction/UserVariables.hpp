/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "CCZ4UserVariables.hpp"

// assign an enum to each variable
enum
{
    c_chi_plot,

    c_Weyl4_Re,
    c_Weyl4_Im,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {"chi",

                                                         "Weyl4_Re",
                                                         "Weyl4_Im"};
}

#endif /* USERVARIABLES_HPP */
