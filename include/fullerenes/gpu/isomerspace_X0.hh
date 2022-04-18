#pragma once
#include "isomerspace_kernel.hh"

class IsomerspaceX0 : public IsomerspaceKernel
{
public:
    void zero_order_geometry(const device_real_t scalerad = 4.0);
    void check_batch();
    
    IsomerspaceX0(const size_t N);
    ~IsomerspaceX0();
};
