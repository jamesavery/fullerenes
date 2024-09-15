#pragma once
#include <CL/sycl.hpp>
#include <numeric>
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <array>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#define FLOAT_TYPEDEFS(T) static_assert(std::is_floating_point<T>::value, "T must be float"); typedef std::array<T,3> coord3d; typedef std::array<T,2> coord2d; typedef T real_t;
#define INT_TYPEDEFS(K) static_assert(std::is_integral<K>::value, "K must be integral type"); typedef std::array<K,3> node3; typedef std::array<K,2> node2; typedef K node_t; typedef std::array<K,6> node6;
#define TEMPLATE_TYPEDEFS(T,K) FLOAT_TYPEDEFS(T) INT_TYPEDEFS(K)

using namespace sycl;
enum ForcefieldType : uint8_t {WIRZ, PEDERSEN, FLATNESS_ENABLED, FLAT_BOND, BOND, ANGLE, DIH, ANGLE_M, ANGLE_P, DIH_A, DIH_M, DIH_P};

enum class EigensolveMode 
{   
    FULL_SPECTRUM,          //Compute all eigenvalues
    ENDS,                   //Compute the smallest and largest non-zero eigenvalues
    FULL_SPECTRUM_VECTORS,  //Compute all eigenvalues and eigenvectors
    ENDS_VECTORS,            //Compute the smallest and largest (non-zero) eigenvalues and corresponding eigenvectors
    SPECIAL                 //Compute the 6 eigenvectors corresponding to the 6 degrees of freedom (3 translations and 3 rotations)
};

enum class LaunchPolicy
{
    ASYNC,
    SYNC
};

enum class IsomerStatus
{
    EMPTY,
    CONVERGED,
    PLZ_CHECK,
    FAILED,
    NOT_CONVERGED
};
template <typename T>
sycl::exception copy(sycl::buffer<T, 1> &dst, sycl::buffer<T, 1> &src)
{
    try
    {
        sycl::host_accessor dst_acc(dst, sycl::write_only);
        sycl::host_accessor src_acc(src, sycl::read_only);
        for (size_t i = 0; i < src.size(); i++)
        {
            dst_acc[i] = src_acc[i];
        }
    }
    catch (sycl::exception e)
    {
        return e;
    }
    return sycl::exception(std::error_code());
}

template <typename T>
sycl::exception copy(sycl::buffer<T, 1> &dst, T *src)
{
    try
    {
        sycl::host_accessor dst_acc(dst, sycl::write_only);
        for (size_t i = 0; i < dst.size(); i++)
        {
            dst_acc[i] = src[i];
        }
    }
    catch (sycl::exception e)
    {
        return e;
    }
    return sycl::exception(std::error_code());
}

template <typename T>
sycl::exception copy(T *dst, sycl::buffer<T, 1> &src)
{
    try
    {
        sycl::host_accessor src_acc(src, sycl::read_only);
        for (size_t i = 0; i < src.size(); i++)
        {
            dst[i] = src_acc[i];
        }
    }
    catch (sycl::exception e)
    {
        return e;
    }
    return sycl::exception(std::error_code());
}
template <typename T, typename K>
struct IsomerBatch
{
    TEMPLATE_TYPEDEFS(T, K);
    private:
        size_t m_capacity = 0;
        size_t m_size = 0;
        size_t n_atoms = 0;
        size_t n_faces = 0;
    public:
        buffer<coord3d, 1> X;
        buffer<coord2d, 1> xys;
        buffer<K, 1> cubic_neighbours;
        buffer<K, 1> dual_neighbours;
        buffer<K, 1> face_degrees;
        buffer<size_t, 1> IDs;
        buffer<size_t, 1> iterations;
        buffer<IsomerStatus, 1> statuses;

        bool allocated = false;

        // std::vector<std::tuple<void**,size_t,bool>> pointers;

        IsomerBatch(size_t n_atoms, size_t n_isomers) : n_atoms(n_atoms), m_capacity(n_isomers),
                                                        n_faces(n_atoms / 2 + 2),
                                                        X                 (range<1>(n_isomers * n_atoms)), 
                                                        xys               (range<1>(n_isomers * n_atoms)), 
                                                        cubic_neighbours  (range<1>(n_isomers * n_atoms * 3)), 
                                                        dual_neighbours   (range<1>(6 * n_isomers * (n_atoms / 2 + 2))), 
                                                        face_degrees      (range<1>((n_atoms / 2 + 2) * n_isomers)), 
                                                        IDs               (range<1>(n_isomers)), 
                                                        iterations        (range<1>(n_isomers)), 
                                                        statuses          (range<1>(n_isomers))
        {   
            sycl::host_accessor X_acc(X, no_init);
            sycl::host_accessor xys_acc(xys, no_init);
            sycl::host_accessor cubic_neighbours_acc(cubic_neighbours, no_init);
            sycl::host_accessor dual_neighbours_acc(dual_neighbours, no_init);
            sycl::host_accessor face_degrees_acc(face_degrees, no_init);
            sycl::host_accessor IDs_acc(IDs, no_init);
            sycl::host_accessor iterations_acc(iterations, no_init);
            sycl::host_accessor statuses_acc(statuses, no_init);




            for (size_t i = 0; i < n_isomers; i++)
            {
                for (size_t j = 0; j < n_atoms; j++)
                {
                      X_acc[i * n_atoms + j] = coord3d{0.0, 0.0, 0.0};
                      xys_acc[i * n_atoms + j] = coord2d{0.0, 0.0};
                    cubic_neighbours_acc[i * n_atoms * 3 + j * 3 + 0] = std::numeric_limits<K>::max();
                    cubic_neighbours_acc[i * n_atoms * 3 + j * 3 + 1] = std::numeric_limits<K>::max();
                    cubic_neighbours_acc[i * n_atoms * 3 + j * 3 + 2] = std::numeric_limits<K>::max();
                }
                for (size_t j = 0; j < 6 * n_faces; j++)
                {
                    dual_neighbours_acc[i * 6 * n_faces + j] = std::numeric_limits<K>::max();
                }
                for (size_t j = 0; j < n_faces; j++)
                {
                    face_degrees_acc[i * n_faces + j] = std::numeric_limits<K>::max();
                }
                IDs_acc[i] = std::numeric_limits<size_t>::max();
                iterations_acc[i] = 0;
                statuses_acc[i] = IsomerStatus::EMPTY;
            }
        }

        size_t size() const {return m_size;}
        size_t capacity() const {return m_capacity;}
        size_t N() const {return n_atoms;}
        size_t Nf() const {return n_faces;}

        IsomerBatch& operator=(const IsomerBatch& other) = default;
        void resize(size_t new_capacity){
            assert(new_capacity >= m_size);
            m_capacity = new_capacity;
            auto new_X = buffer<coord3d, 1>( range<1>(m_capacity * n_atoms)); copy(new_X, X); X = new_X;
            auto new_xys = buffer<coord2d, 1>( range<1>(m_capacity * n_atoms)); copy(new_xys, xys); xys = new_xys;
            auto new_cubic_neighbours = buffer<K, 1>( range<1>(m_capacity * n_atoms * 3)); copy(new_cubic_neighbours, cubic_neighbours); cubic_neighbours = new_cubic_neighbours;
            auto new_dual_neighbours = buffer<K, 1>( range<1>(6 * m_capacity * n_faces)); copy(new_dual_neighbours, dual_neighbours); dual_neighbours = new_dual_neighbours;
            auto new_face_degrees = buffer<K, 1>( range<1>(m_capacity * n_faces)); copy(new_face_degrees, face_degrees); face_degrees = new_face_degrees;
            auto new_IDs = buffer<size_t, 1>( range<1>(m_capacity)); copy(new_IDs, IDs); IDs = new_IDs;
            auto new_iterations = buffer<size_t, 1>( range<1>(m_capacity)); copy(new_iterations, iterations); iterations = new_iterations;
            auto new_statuses = buffer<IsomerStatus, 1>( range<1>(m_capacity)); copy(new_statuses, statuses); statuses = new_statuses;
        }
        ~IsomerBatch() = default;
};



template<typename T, typename K>
void copy(IsomerBatch<T,K>& dst, const IsomerBatch<T,K>& src){
    assert(dst.N() == src.N());
    assert(dst.capacity() == src.capacity());
    copy(dst.cubic_neighbours, src.cubic_neighbours);
    copy(dst.dual_neighbours, src.dual_neighbours);
    copy(dst.statuses, src.statuses);
    copy(dst.X, src.X);
    copy(dst.xys, src.xys);
    copy(dst.iterations, src.iterations);
    copy(dst.IDs, src.IDs);
}

template <typename T, typename K>
void fill(FullereneBatch<T,K>& B, int mytask_id = 0, int ntasks = 1) {
  int N = B.N_;
  int Nf = B.Nf_;
  int N_graphs = B.capacity();
  auto face_degrees_acc = B.d_.deg_;
  auto dual_neighbours_acc = B.d_.A_cubic_;
  auto statuses_acc = B.m_.flags_;
  BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,false, false, mytask_id, ntasks);
  Graph G;

  G.neighbours = neighbours_t(Nf, std::vector<node_t>(6,-1));
  G.N = Nf;
  int num_generated = 0;
  for (int i = 0; i < N_graphs; ++i) {
    bool more_isomers = BuckyGen::next_fullerene(BuckyQ, G);
    if (!more_isomers) break;
    num_generated++;
    statuses_acc[i] = StatusFlag::DUAL_INITIALIZED;
    for(int j = 0; j < Nf; j++) {
      face_degrees_acc[i*Nf + j] = G.neighbours[j].size();
      for(int k = 0; k < G.neighbours[j].size(); k++) 
        dual_neighbours_acc[i*Nf*6 + j*6 + k] = G.neighbours[j][k];

    }
  }
  BuckyGen::stop(BuckyQ);
  if (num_generated < N_graphs) 
    for (int i = num_generated; i < N_graphs; ++i) {
      statuses_acc[i] = StatusFlag::DUAL_INITIALIZED;
      //Repeat the same graphs as already generated.
      for(int j = 0; j < Nf; j++) {
        face_degrees_acc[i*Nf + j] = face_degrees_acc[(i%num_generated)*Nf + j];
        for(int k = 0; k < 6; k++) 
          dual_neighbours_acc[i*Nf*6 + j*6 + k] = dual_neighbours_acc[(i%num_generated)*Nf*6 + j*6 + k];

      }
    }
}