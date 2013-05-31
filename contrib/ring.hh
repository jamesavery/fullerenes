#ifndef RING_HH
#define RING_HH
#include <inttypes.h>

template <typename T>
class ring {
    public:

        virtual ~ring() {}
        virtual T add(T a, T b) = 0;
        virtual T sub(T a, T b) = 0;
        virtual T mul(T a, T b) = 0;
        virtual T add_inv(T a) = 0;
        virtual T add_id() = 0;
        virtual T mul_id() = 0;
};

class ringZ : public ring<int64_t> {
    public:
        virtual int64_t add(int64_t a, int64_t b) { return a + b; }
        virtual int64_t sub(int64_t a, int64_t b) { return a - b; }
        virtual int64_t mul(int64_t a, int64_t b) { return a * b; }
        virtual int64_t add_inv(int64_t a)    { return -a; }
        virtual int64_t add_id()          { return 0; }
        virtual int64_t mul_id()          { return 1; }
};

#endif
