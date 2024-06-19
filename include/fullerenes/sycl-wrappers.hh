#include <string>
#include <array>
#include "graph.hh"
#include "polyhedron.hh"
/* enum class IsomerStatus
{
    EMPTY,
    CONVERGED,
    PLZ_CHECK,
    FAILED,
    NOT_CONVERGED
}; */

enum class FullereneGraphType
{
    CUBIC,
    DUAL,
    BOTH
};

enum class DeviceType{
    CPU,
    GPU,
    ACCELERATOR,
    HOST
};

enum class DeviceProperty{
    MAX_WORK_GROUP_SIZE,
    MAX_CLOCK_FREQUENCY,
    MAX_COMPUTE_UNITS,
    MAX_MEM_ALLOC_SIZE,
    GLOBAL_MEM_SIZE,
    LOCAL_MEM_SIZE,
    MAX_CONSTANT_BUFFER_SIZE,
    MAX_CONSTANT_ARGS,
    NUMBER_OF_PROPERTIES
};

template <typename T>
struct SyclVector
{
    SyclVector();
    SyclVector(size_t size);
    SyclVector(size_t size, T value);
    SyclVector(const SyclVector<T>& other);
    SyclVector(SyclVector<T>&& other);
    SyclVector<T>& operator=(const SyclVector<T>& other);

    ~SyclVector();

    void fill(T data);
    T* data();
    size_t size();
    size_t capacity();
    void resize(size_t new_size);
    void reserve(size_t new_capacity);
    void clear();

    T& operator[](size_t index);
    const T& operator[](size_t index) const;

    T& at(size_t index);
    const T& at(size_t index) const;

    void push_back(const T& value);
    void pop_back();
    
    private:
        size_t size_;
        size_t capacity_;
        T* data_;
};

template <typename T = float, typename K = uint16_t>
struct FullereneIsomer{
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    static_assert(std::is_integral<K>::value, "K must be an integral type");

    FullereneIsomer();
    ~FullereneIsomer();

    FullereneIsomer(size_t N, 
                    bool allocate_hessian_etc = false);
    
    FullereneIsomer(const Graph& G, bool is_cubic = false, bool allocate_hessian_etc = false);
    FullereneIsomer(const neighbours_t& neighbours, bool is_cubic = false, bool allocate_hessian_etc = false);
    FullereneIsomer(const Polyhedron& P, bool is_cubic = true, bool allocate_hessian_etc = false);
    FullereneIsomer(const FullereneIsomer<T,K>& other);
    FullereneIsomer(FullereneIsomer<T,K>&& other);
    FullereneIsomer<T,K>& operator=(const FullereneIsomer<T,K>& other);


    SyclVector<T> X_cubic_;      //3D Embedding of the Cubic Graph {N * 3}
    SyclVector<T> X_dual_;       //3D Embedding of the Dual Graph {Nf * 6}
    SyclVector<K> A_cubic_;      //Adjacency Matrix (Cubic) {N * 3}
    SyclVector<K> A_dual_;       //Adjacency Matrix (Dual) {Nf * 6}
    SyclVector<K> faces_;        //Atom indices of the faces {Nf * 6}
    SyclVector<K> deg_;          //Vertex degrees in the dual graph, face degrees in the cubic graph, face degrees in the dual graph is always 3 (it is a triangulation)
    SyclVector<T> gradient_;     //Forcefield gradient {N * 3}

    //Not Allocated by default
    SyclVector<T> hessian_;      //Forcefield hessian {N * 90}
    SyclVector<T> eigenvalues_;  //Eigenvalues of the hessian matrix {N * 3}
    SyclVector<T> eigenvectors_; //Eigenvectors of the hessian matrix {(N*3)^2}

    size_t N_;                   //Number of vertices in the cubic graph {1}
    size_t Nf_;                  //Number of faces in the dual graph {1}
    //IsomerStatus status;        //Convergence status of the isomer {1}
    size_t ID_;                  //Buckygen ID of the isomer {1}           
    size_t iterations_;          //Number of forcefield CG iterations {1}
    T energy_;                   //Forcefield energy of the isomer {1}
    T eccentricity_;             //Eccentricity of the isomer {1}
    T volume_divergence_;        //Volume divergence of the isomer {1}
    FullereneGraphType g_type_;  //Graph Representation, Cubic, Dual or Both {1}
    std::array<K,12> spiral_;    //Spiral code, named canonically after pentagon indices {1*12}
};



struct DeviceWrapper
{   
    DeviceWrapper();
    DeviceWrapper(size_t device_id, const DeviceType device_type);
    ~DeviceWrapper() = default;
    bool is_cpu() const;
    bool is_gpu() const;
    size_t get_id() const;
    std::string get_name() const;
    size_t get_property(DeviceProperty property) const;

    size_t id_;
    DeviceType type_;
    const void* device_;

    static std::vector<std::byte> __FULLERENE_ALLCPUDEVICES__;
    static std::vector<std::byte> __FULLERENE_ALLGPUDEVICES__;
    static std::vector<std::byte> __FULLERENE_ALLACCELERATORS__;
};

struct QueueWrapper
{
    QueueWrapper(const int device_id = 0, const DeviceType = DeviceType::GPU, bool in_order = true);
    ~QueueWrapper() = default;
    void wait();
    const DeviceWrapper get_device() const;

    size_t id_;
    const DeviceWrapper device_;
    std::vector<std::byte> queue_;
    //static std::vector<std::byte> allocator_queue_;
};