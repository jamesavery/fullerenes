#pragma once
// Single-definition shim. The SyclEventImpl / SyclQueueImpl structs and the
// SyclQueue / SyclEvent / Device method bodies live ONCE in queue-impl.hh;
// defining DEFINE_SYCL_QUEUE_METHODS before including it emits the method
// definitions in this translation unit. (Included by kernel.cc; the production
// library compiles the same methods via sycl-util-impl.cc, which uses this idiom.)
#define DEFINE_SYCL_QUEUE_METHODS
#include "queue-impl.hh"
