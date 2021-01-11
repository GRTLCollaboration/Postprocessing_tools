/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef REPROCESSINGLEVEL_HPP_
#define REPROCESSINGLEVEL_HPP_

#include "GRAMRLevel.hpp"
#ifdef USE_AHFINDER
#include "AHInterpolation.hpp"
#include "AHSphericalCoords.hpp"
#include "ApparentHorizon.hpp"
#endif

class ReprocessingLevel : public GRAMRLevel {
  friend class DefaultLevelFactory<ReprocessingLevel>;
  // Inherit the contructors from GRAMRLevel
  using GRAMRLevel::GRAMRLevel;
  BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);

  // initialize data
  virtual void initialData() { m_state_new.setVal(0.); }

  void postRestart() {
    // Add code here to do what you need it to do on each level
    pout() << "The time is " << m_time << " on level " << m_level
           << ". Your wish is my command." << endl;

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
      m_bh_amr.ah_finder.solve(m_dt, m_time, m_restart_time);
#endif
  }

  virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                               const double a_time) {}

  virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                       const FArrayBox &current_state) {
    tagging_criterion.setVal(0.);
  };
};

#endif /* REPROCESSINGLEVEL_HPP_ */
