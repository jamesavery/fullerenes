#pragma once
#include "isomerspace_kernel.hh"
#include "gpudatastruct.hh"

class IsomerspaceTutte : public IsomerspaceKernel<FullereneGraph>
{
public:
    struct IsomerBatchStats : GPUDataStruct{
        size_t*       iteration_counts;
        IsomerStatus* isomer_statuses;
        size_t*       isomer_IDs;
        uint8_t*      Nface; //Face type of the outer face (for fullerenes Either 5 or 6);

        IsomerBatchStats(){
            pointers = {{"iteration_counts", (void**)&iteration_counts, sizeof(size_t)}, 
                        {"isomer_statuses", (void**)&isomer_statuses, sizeof(IsomerStatus)}, 
                        {"isomer_IDs", (void**)&isomer_IDs, sizeof(size_t)},
                        {"Nface", (void**)&Nface, sizeof(uint8_t)}};
        }
    };

    struct IsomerBatch : GPUDataStruct{
        device_real_t* xys;        //Nx2xBatchSize     
        device_node_t* outer_face; //(5-6)xBatchSize : Size of this array is variable since the individual vectors are of unknown size (5 or 6 for fullerenes though)
        device_node_t* neighbours; //Nx3xBatchSize

        IsomerBatchStats stats;

        IsomerBatch(){
            pointers =  {{"xys",(void**)&xys,sizeof(device_real_t)*2}, {"neighbours",(void**)&neighbours, sizeof(device_node_t)*3}, 
                        {"outer_face", (void**)&outer_face, sizeof(device_node_t)*1}};
        }
    };

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

    std::vector<IsomerBatch> d_batch;
    std::vector<IsomerBatch> h_batch;

};


