/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef REPROCESSINGLEVEL_HPP_
#define REPROCESSINGLEVEL_HPP_

#include "GRAMRLevel.hpp"
#include "WeylExtraction.hpp"

class ReprocessingLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ReprocessingLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData() { m_state_new.setVal(0.); }

    void postRestart()
    {
        // Add code here to do what you need it to do on each level
        pout() << "The time is " << m_time << " on level " << m_level
               << ". Your wish is my command." << endl;

        // don't forget to setup the interpolator in ReprocessingTool.cpp
        // and don't add a 'restart_file' parameter so that 'setupAMRObject'
        // doesn't restart before the interpolator is set
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            // Now refresh the interpolator and do the interpolation
            fillAllGhosts();
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                         m_time == 0., m_restart_time);
            my_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state)
    {
        tagging_criterion.setVal(0.);
    };
};

#endif /* REPROCESSINGLEVEL_HPP_ */
