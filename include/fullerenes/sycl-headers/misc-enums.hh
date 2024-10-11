#include <array>
#include <cstdint>

#define FLOAT_TYPEDEFS(T) static_assert(std::is_floating_point<T>::value, "T must be float"); typedef std::array<T,3> coord3d; typedef std::array<T,2> coord2d; typedef T real_t;
#define INT_TYPEDEFS(K) static_assert(std::is_integral<K>::value, "K must be integral type"); typedef std::array<K,3> node3; typedef std::array<K,2> node2; typedef K node_t; typedef std::array<K,6> node6;
#define TEMPLATE_TYPEDEFS(T,K) FLOAT_TYPEDEFS(T) INT_TYPEDEFS(K)

using uint8_t = unsigned char;
using uint16_t = unsigned short;
using uint32_t = unsigned int;
using uint64_t = unsigned long;

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