#include <fullerenes/sycl-isomer-batch.hh>

int main() {
    IsomerBatch<float,uint16_t> batch(20, 100);
    IsomerBatch<float,uint8_t> batch2(20, 100);
    IsomerBatch<float,uint32_t> batch3(20, 100);
    IsomerBatch<float,uint64_t> batch4(20, 100);
    IsomerBatch<float,size_t> batch5(20, 100);
    IsomerBatch<float,int> batch6(20, 100);
    IsomerBatch<float,unsigned int> batch7(20, 100);
    IsomerBatch<float,long> batch8(20, 100);
    IsomerBatch<float,unsigned long> batch9(20, 100);
    IsomerBatch<float,long long> batch10(20, 100);
    IsomerBatch<float,unsigned long long> batch11(20, 100);
    IsomerBatch<float,short> batch12(20, 100);
    IsomerBatch<float,unsigned short> batch13(20, 100);
    
    IsomerBatch<double,uint16_t> batch14(20, 100);
    IsomerBatch<double,uint8_t> batch15(20, 100);
    IsomerBatch<double,uint32_t> batch16(20, 100);
    IsomerBatch<double,uint64_t> batch17(20, 100);
    IsomerBatch<double,size_t> batch18(20, 100);
    IsomerBatch<double,int> batch19(20, 100);
    IsomerBatch<double,unsigned int> batch20(20, 100);
    IsomerBatch<double,long> batch21(20, 100);
    IsomerBatch<double,unsigned long> batch22(20, 100);
    IsomerBatch<double,long long> batch23(20, 100);
    IsomerBatch<double,unsigned long long> batch24(20, 100);
    IsomerBatch<double,short> batch25(20, 100);
    IsomerBatch<double,unsigned short> batch26(20, 100);
    
    return 0;
}