/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "EmptyDiagnosticVariables.hpp"

// assign an enum to each variable
enum
{
    c_chi,

    c_Ham,

    c_Weyl4_Re,
    c_Weyl4_Im,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {
    "chi", "Ham", "Weyl4_Re", "Weyl4_Im"};
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
