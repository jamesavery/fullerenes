#pragma once
#include "isomerspace_kernel.hh"

class IsomerspaceTutte : public IsomerspaceKernel
{
public:

    void tutte_layout();
    void check_batch();

    IsomerspaceTutte(const size_t N);
    ~IsomerspaceTutte();

protected:
    const size_t TUTTE_MAX_ITERATION = 1000;
};

