#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/argparser.hh>
#include <fullerenes/linalg/qr/batched_qr.hh>


int main(int argc, char** argv){
    //Tall and skinny matrices
    SyclQueue ctx(std::string("gpu"));
    int m = argc > 1 ? std::stoi(argv[1]) : 100;
    int n = argc > 2 ? std::stoi(argv[2]) : 10;
    int batch_size = argc > 3 ? std::stoi(argv[3]) : ctx.device().get_property(DeviceProperty::MAX_COMPUTE_UNITS);

    auto A = SyclVector<float>(m*n*batch_size);
    auto workspace = SyclVector<std::byte>(chol_qr_batched_buffer_size<float>(ctx, false, batch_size, m*n, m, n, A.to_span()));

    std::cout << "Workspace made" << std::endl;

    std::transform(A.begin(), A.end(), A.begin(), [](auto val){return (float) (rand()%100);});

    std::cout << "Matrix made" << std::endl;

    auto T0 = std::chrono::high_resolution_clock::now();
    chol_qr_batched<float>( ctx, 
                            false, 
                            batch_size, 
                            m*n, 
                            m, 
                            n, 
                            A.to_span(), 
                            workspace.to_span());
    auto T1 = std::chrono::high_resolution_clock::now();


    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(T1 - T0).count()/batch_size << " µs / fullerene" << std::endl;
}