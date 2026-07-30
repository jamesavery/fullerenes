#pragma once

// =====================================================================
// optim::SeqCtx -- the sequential execution context.
//
// Every framework entry point takes a Ctx as its first argument (the
// parallel-primitives idiom), so the interface is already shaped for
// the par::Seq / par::Omp / device lowering of migration step 6.  This
// landing is CPU-sequential only: bodies are plain loops and SeqCtx is
// an empty tag threaded through signatures.  Lowering the bodies onto
// the primitives (with D1 fixed-shape reductions) is its own validated
// step, not a recompile -- see PRIMITIVES-SPEC.md 10.
// =====================================================================

namespace optim {

struct SeqCtx {};

}  // namespace optim
