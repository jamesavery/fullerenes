#pragma once
#include "isomerspace_kernel.hh"
#include "gpudatastruct.hh"

class IsomerspaceTutte : public IsomerspaceKernel<FullereneGraph>
{
public:

    void tutte_layout();
    void update_batch();
    void check_batch();
    void eject_isomer(size_t i, size_t idx);

    IsomerspaceTutte(const size_t N);
    ~IsomerspaceTutte();

protected:
    size_t converged_count = 0;             //Total number of converged isomers optimized by this object.
    size_t failed_count = 0;                //Total number of failed isomers optimized by this object.  
    const size_t TUTTE_MAX_ITERATION = 1000;
};


