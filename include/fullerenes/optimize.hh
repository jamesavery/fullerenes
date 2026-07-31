#pragma once

// =====================================================================
// fullerenes/optimize.hh -- the unified optimizer framework (umbrella).
//
// Graduated from claude-projects/optimize (first landing 2026-07-19,
// graduation 2026-07-30); that sub-project remains the validation
// harness (re-expression gates, FD pins, head-to-head validators).
//
// A globalization PARADIGM (LineSearch / TrustRegion) owns the outer
// iteration loop; the STEP (LBFGS, SteihaugCG, LevenbergBisect, ...) is
// a policy valid only within its paradigm (compile-time-checked); the
// EnergyModel advertises its faces through capability concepts.  See
// DESIGN.md and the per-header comments.
//
// This umbrella documents the FULL public surface and therefore drags
// in every model and driver (Deltahedron machinery via delta_optimize,
// the Alexandrov DCEL via alex_polish).  Consumers wanting one path
// should include its per-path headers instead -- e.g. a line-search
// L-BFGS user needs only optim/linesearch.hh + optim/steps/lbfgs.hh
// (+ their model header).
// =====================================================================

#include "fullerenes/optim/core.hh"
#include "fullerenes/optim/seq_ctx.hh"
#include "fullerenes/optim/linalg.hh"
#include "fullerenes/optim/preconditioner.hh"

#include "fullerenes/optim/steps/directions.hh"
#include "fullerenes/optim/steps/lbfgs.hh"
#include "fullerenes/optim/steps/levenberg.hh"
#include "fullerenes/optim/steps/steihaug.hh"

#include "fullerenes/optim/linesearch.hh"
#include "fullerenes/optim/trustregion.hh"
#include "fullerenes/optim/driver.hh"

#include "fullerenes/optim/models/extwu_angle.hh"
#include "fullerenes/optim/models/geometry_hessians.hh"
#include "fullerenes/optim/models/alex_kappa.hh"
#include "fullerenes/optim/models/aet.hh"

#include "fullerenes/optim/alex_polish.hh"
#include "fullerenes/optim/delta_optimize.hh"
