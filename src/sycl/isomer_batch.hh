#include "CL/sycl.hpp"

using namespace sycl;

#define UINT_TYPE uint16_t
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
                                                      face_degrees      (range<1>((n_atoms / 2 + 2) * 1)), 
                                                      IDs               (range<1>(n_isomers)), 
                                                      iterations        (range<1>(n_isomers)), 
                                                      statuses          (range<1>(n_isomers))
      {   
          X                 = buffer<coord3d, 1>( range<1>(n_isomers * n_atoms));
          xys               = buffer<coord2d, 1>( range<1>(n_isomers * n_atoms));
          cubic_neighbours  = buffer<K, 1>( range<1>(n_isomers * n_atoms * 3));
          dual_neighbours   = buffer<K, 1>( range<1>(6 * n_isomers * n_faces));
          face_degrees      = buffer<K, 1>( range<1>(n_isomers * n_faces));
          IDs               = buffer<size_t, 1>( range<1>(n_isomers));
          iterations        = buffer<size_t, 1>( range<1>(n_isomers));
          statuses          = buffer<IsomerStatus, 1>( range<1>(n_isomers));
          
          host_accessor X_acc(X, no_init);
          host_accessor xys_acc(xys, no_init);
          host_accessor cubic_neighbours_acc(cubic_neighbours, no_init);
          host_accessor dual_neighbours_acc(dual_neighbours, no_init);
          host_accessor face_degrees_acc(face_degrees, no_init);
          host_accessor IDs_acc(IDs, no_init);
          host_accessor iterations_acc(iterations, no_init);
          host_accessor statuses_acc(statuses, no_init);

          


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


};

