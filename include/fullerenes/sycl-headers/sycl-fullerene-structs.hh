#pragma once
#include <fullerenes/sycl-headers/sycl-status-enum.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-data.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-misc-tuple-fun.hh>
#include <fullerenes/sycl-headers/sycl-fullerene.hh>
#include <fullerenes/polyhedron.hh>
#include <bitset>
#include <bit>
#include <thread>
#include <ranges>
#include <source_location>
using uint16_t = unsigned short;

#define OOR_ERROR(index, max) std::out_of_range(std::string(__FILE__) + ":" + std::to_string(__LINE__) + " Index: " + std::to_string(index) + " out of range: 0 - " + std::to_string(max))

template <typename T, typename K>
struct FullereneBatch;

template <typename T, typename K>
struct FullereneQueue;

template <typename T, typename K>
struct FullereneBatchView;

template <typename T, typename K>
struct Fullerene;

//Encapsulation of functions that must be friends of FullereneBatch and FullereneQueue so as to modify, size_, capacity_, front_, back_.
struct QueueUtil{
    template<typename T1, typename K1>
    static void push(SyclQueue& Q, FullereneQueue<T1, K1> &dst_queue, FullereneBatch <T1, K1>&src_batch, ConditionFunctor transfer_cond = ConditionFunctor(), StatusEnum consumed_status = StatusEnum(0));
    
    template<typename T1, typename K1>
    static void push(SyclQueue& Q, FullereneBatch<T1, K1> &dst_batch, FullereneQueue<T1, K1> &src_queue, ConditionFunctor transfer_cond = ConditionFunctor(), StatusEnum consumed_status = StatusEnum(0));

    template <typename DstBatch, typename SrcBatch>
    static void push_impl(SyclQueue& Q, DstBatch& dst_batch, SrcBatch& src_batch, ConditionFunctor condition = ConditionFunctor(), StatusEnum consumed_status = StatusEnum(0));
};


template <typename T = float, typename K = uint16_t>
struct FullereneBatch
{
    using value_type = Fullerene<T, K>;

    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    static_assert(std::is_integral<K>::value, "K must be an integral type");

    FullereneBatch();
    ~FullereneBatch() = default;

    FullereneBatch(size_t N, int capacity);

    void push_back(const Graph &G, const int ID = -1);
    void push_back(const neighbours_t &neighbours, const int ID = -1);
    void push_back(const PlanarGraph &G, const int ID = -1);
    void push_back(const Polyhedron &P, const int ID = -1);
    void push_back(const Fullerene<T, K> &fullerene, const int ID = -1);
    void resize(size_t new_capacity);

    FullereneBatch(const FullereneBatch<T, K> &other) = default;
    FullereneBatch(FullereneBatch<T, K> &&other) = default;
    FullereneBatch<T, K> &operator=(const FullereneBatch<T, K> &other) = default;

    Fullerene<T, K> operator[](size_t index) const;
    bool operator==(const FullereneBatch<T, K> &other) const {return d_ == other.d_ && m_ == other.m_;}
    constexpr inline operator FullereneBatchView<T, K>() const { return FullereneBatchView<T, K>(*this); }

    template<typename T1, typename K1>
    friend void QueueUtil::push(SyclQueue& Q, FullereneQueue<T1, K1> &dst_queue, FullereneBatch <T1, K1>&src_batch, ConditionFunctor transfer_cond, StatusEnum consumed_status);
    
    template<typename T1, typename K1>
    friend void QueueUtil::push(SyclQueue& Q, FullereneBatch<T1, K1> &dst_batch, FullereneQueue<T1, K1> &src_queue, ConditionFunctor transfer_cond, StatusEnum consumed_status);

    template <typename DstBatch, typename SrcBatch>
    friend void QueueUtil::push_impl(SyclQueue& Q, DstBatch& dst_batch, SrcBatch& src_batch, ConditionFunctor condition, StatusEnum consumed_status);

    // Getters for private members that are not supposed to be modified outside the class
    constexpr inline int size()     const {return size_;}
    constexpr inline int capacity() const {return capacity_;}
    constexpr inline Fullerene<T, K> front()    const {return (*this)[0];}
    constexpr inline Fullerene<T, K> back()     const {return (*this)[size_-1];}

    constexpr inline auto begin() const { return FullereneBatchView<T,K>(*this).begin(); }
    constexpr inline auto end() const { return FullereneBatchView<T,K>(*this).end(); }

    constexpr inline void clear() { size_ = 0; }
    FullereneDataMembers<SyclVector, T, K> d_;
    FullereneMetaMembers<SyclVector, K> m_;

    size_t N_;                      // Number of vertices in the cubic graph {1}
    size_t Nf_;                     // Number of faces in the dual graph {1}

    //Output stream operator for FullereneBatch
    template <typename U, typename V>
    friend std::ostream& operator<<(std::ostream& os, const FullereneBatch<U, V>& batch);

protected:
    int capacity_ = 0; // Maximum number of isomers in the batch {1}
    int size_ = 0;     // Number of isomers in the batch {1}

    void prepare_for_push_back(const neighbours_t &neighbours);
};

template <typename T = float, typename K = uint16_t>
struct QueueIterator
{
    using difference_type = std::ptrdiff_t;
        using value_type =  Fullerene<T,K>; 
        using pointer =  Fullerene<T,K>*;
        using reference =  Fullerene<T,K>&;
    using iterator_category = std::random_access_iterator_tag;

    //Converts the circular index to the index in the batch
    int batch_ix_(int index) const {  return (queue_.get().front_index() + index) % queue_.get().capacity(); } 

    QueueIterator<T,K>(const ReferenceWrapper<const FullereneQueue<T,K>>& queue, int index) : index_(index), queue_(queue) {
        if (index < 0 || index > queue.get().size()) {throw OOR_ERROR(index, queue.get().size());}
    }
    
    QueueIterator<T,K>(const QueueIterator<T,K> &other) = default;
    QueueIterator<T,K> &operator=(const QueueIterator<T,K> &other) = default;

    constexpr inline Fullerene<T,K> operator*() const { return queue_.get()[batch_ix_(index_)]; }
    constexpr inline Fullerene<T,K> operator->() const { return queue_.get()[batch_ix_(index_)]; }
    constexpr inline Fullerene<T,K> operator[](int n) const { return queue_.get()[batch_ix_(index_ + n)]; }

    difference_type operator-(const QueueIterator<T,K> &other) const { return index_ - other.index_; }

    QueueIterator<T,K> &operator++() { index_++; return *this; }
    QueueIterator<T,K> operator++(int) { return QueueIterator<T,K>(queue_, index_++); }
    QueueIterator<T,K> &operator--() { index_--; return *this; }
    QueueIterator<T,K> operator--(int) { return QueueIterator<T,K>(queue_, index_--); }
    QueueIterator<T,K> &operator+=(int n) { index_ += n; return *this; }
    QueueIterator<T,K> &operator-=(int n) { index_ -= n; return *this; }
    QueueIterator<T,K> operator+(int n) const { return QueueIterator<T,K>(queue_, index_ + n); }
    QueueIterator<T,K> operator-(int n) const { return QueueIterator<T,K>(queue_, index_ - n); }

    friend bool operator<(const QueueIterator<T,K> &lhs, const QueueIterator<T,K> &rhs)  { return lhs.index_ < rhs .index_; }
    friend bool operator>(const QueueIterator<T,K> &lhs, const QueueIterator<T,K> &rhs)  { return lhs.index_ > rhs.index_; }
    friend bool operator<=(const QueueIterator<T,K> &lhs, const QueueIterator<T,K> &rhs) { return lhs.index_ <= rhs.index_; }
    friend bool operator>=(const QueueIterator<T,K> &lhs, const QueueIterator<T,K> &rhs) { return lhs.index_ >= rhs.index_; }
    friend bool operator==(const QueueIterator<T,K> &lhs, const QueueIterator<T,K> &rhs) { return lhs.index_ == rhs.index_; }
    friend bool operator!=(const QueueIterator<T,K> &lhs, const QueueIterator<T,K> &rhs) { return lhs.index_ != rhs.index_; }


    int get_index() const { return index_; }
private:
    int index_;
    ReferenceWrapper<const FullereneQueue<T,K>> queue_;
};

// FullereneQueue and FullereneBatch are functionally the we use the type system to express intent and to prevent misuse,
// FullereneQueues have the additional property that their non-StatusEnum::EMPTY isomers are contiguously stored in memory
using namespace condition_detail;
template <typename T = float, typename K = uint16_t>
struct FullereneQueue : public FullereneBatch<T,K>
{
    using FullereneBatch<T,K>::d_;
    using FullereneBatch<T,K>::m_;
    using FullereneBatch<T,K>::N_;
    using FullereneBatch<T,K>::Nf_;
    using FullereneBatch<T,K>::capacity_;
    using FullereneBatch<T,K>::size_;


    FullereneQueue() = default;
    ~FullereneQueue() = default;

    FullereneQueue(size_t N, int capacity) : FullereneBatch<T, K>(N, capacity), front_(-1), back_(-1) {}
    FullereneQueue(const FullereneQueue<T, K> &other) = default;
    FullereneQueue(FullereneQueue<T, K> &&other) = default;
    FullereneQueue<T, K>& operator=(const FullereneQueue<T, K> &other) = default;
    FullereneQueue<T, K>& operator=(FullereneQueue<T, K> &&other) = default;

    template<typename T1, typename K1>
    friend void QueueUtil::push(SyclQueue& Q, FullereneQueue<T1, K1> &dst_queue, FullereneBatch <T1, K1>&src_batch, ConditionFunctor transfer_cond, StatusEnum consumed_status);
    
    template<typename T1, typename K1>
    friend void QueueUtil::push(SyclQueue& Q, FullereneBatch<T1, K1> &dst_batch, FullereneQueue<T1, K1> &src_queue, ConditionFunctor transfer_cond, StatusEnum consumed_status);

    template <typename DstBatch, typename SrcBatch>
    friend void QueueUtil::push_impl(SyclQueue& Q, DstBatch& dst_batch, SrcBatch& src_batch, ConditionFunctor condition, StatusEnum consumed_status);

    Fullerene<T, K> operator[](size_t index) const {
        if (index >= this->size()) {throw OOR_ERROR(index, this->size());}
        return (FullereneBatchView<T, K>(*this, 0, this->capacity()))[batch_ix(index)];
    }

    void push_back(const Graph &G, const int ID = -1, const StatusEnum status = StatusEnum(0)) {
        queue_prepare_push_back(G.neighbours);
        (*this)[size_ - 1] = static_cast<std::tuple<std::reference_wrapper<const Graph>, size_t>>(std::tuple{std::ref(G), ID});
        (*this)[size_ - 1].m_.flags_.get() |= status;
    }
    
    void push_back(const Polyhedron& P, const int ID = -1, const StatusEnum status = StatusEnum(0)){
        queue_prepare_push_back(P.neighbours);
        (*this)[size_ - 1] = static_cast<std::tuple<std::reference_wrapper<const Polyhedron>, size_t>>(std::tuple{std::ref(P), ID});
        (*this)[size_ - 1].m_.flags_.get() |= status;
    }
    
    int batch_ix(int index) const { return (front_ + index) % capacity_;}

    constexpr inline auto begin() const { return QueueIterator<T,K>( std::ref(*this), 0) ; }
    constexpr inline auto end() const { return QueueIterator<T,K>( std::ref(*this), this->size()); }
    constexpr inline bool empty() const { return this->size() == 0; }
    constexpr inline auto front() const { if (empty()) {throw std::out_of_range("Queue is empty");} return (*this)[front_]; }
    constexpr inline auto back() const { if (empty()) {throw std::out_of_range("Queue is empty");} return (*this)[back_]; }
    constexpr inline auto front_index() const { return front_; }
    constexpr inline auto back_index() const { return back_; }
    constexpr inline int size() const { return size_; }

    void resize(size_t new_capacity);
    
private:
    int front_ = -1;    // Index of the first isomer in the queue {1}
    int back_ = -1;     // Index of the last isomer in the queue {1}

    void queue_prepare_push_back(const neighbours_t &neighbours) {
        if (N_ == 0 || Nf_ == 0) { 
            auto N = neighbours.size();
            bool is_cubic = neighbours[0].size() == 3;
            N_ = is_cubic ? N : (N - 2) * 2;
            Nf_ = is_cubic ? N/2 + 2 : N;
        }
        size_++;
        if ((size_-1) == capacity_) {
            resize(capacity_ == 0 ? 1 : capacity_ * 2);
        } else if (size_ == 1) {
            front_ = 0;
            back_ = 0;
        } else {
            back_ = (back_ + 1) % capacity_;
        }
    }
};

template <typename T = float, typename K = uint16_t>
struct FullereneBatchView
{
    using value_type = Fullerene<T, K>;
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    static_assert(std::is_integral<K>::value, "K must be an integral type");

    struct Iterator
    {   
        using iterator_category = std::random_access_iterator_tag;
        using difference_type = int;
        using pointer = int;
        using value_type = Fullerene<T,K>;
        using reference = Fullerene<T,K>;


        Iterator(const FullereneBatchView<T,K>& view, int index) : index_(index), view_(view) {
            if (index < 0 || index > view.size()) {throw OOR_ERROR(index, view.size());}
        }
        Iterator(const Iterator &other) = default;
        Iterator(Iterator &&other) = default;
        Iterator &operator=(const Iterator &other) = default;
        Iterator &operator=(Iterator &&other) = default;

        Fullerene<T,K> operator*() const { return view_[index_]; }
        Fullerene<T,K> operator->() const { return view_[index_]; }
        Fullerene <T,K> operator[](int n) const { return view_[index_ + n]; }

        difference_type operator-(const Iterator &other) const { return index_ - other.index_; }

        Iterator &operator++() { index_++; return *this; }
        Iterator operator++(int) { return Iterator(view_, index_++); }
        Iterator &operator--() { index_--; return *this; }
        Iterator operator--(int) { return Iterator(view_, index_--); }
        Iterator &operator+=(int n) { index_ += n; return *this; }
        Iterator &operator-=(int n) { index_ -= n; return *this; }
        Iterator operator+(int n) const { return Iterator(view_, index_ + n); }
        Iterator operator-(int n) const { return Iterator(view_, index_ - n); }


        friend bool operator<(const Iterator &lhs, const Iterator &rhs)  { return lhs.index_ < rhs.index_; }
        friend bool operator>(const Iterator &lhs, const Iterator &rhs)  { return lhs.index_ > rhs.index_; }
        friend bool operator<=(const Iterator &lhs, const Iterator &rhs) { return lhs.index_ <= rhs.index_; }
        friend bool operator>=(const Iterator &lhs, const Iterator &rhs) { return lhs.index_ >= rhs.index_; }
        friend bool operator==(const Iterator &lhs, const Iterator &rhs) { return lhs.index_ == rhs.index_; }
        friend bool operator!=(const Iterator &lhs, const Iterator &rhs) { return lhs.index_ != rhs.index_; }
        
        int get_index() const { return index_; }
        auto get_view() const { return view_; }
        private:
            pointer index_;
            const FullereneBatchView<T,K> view_;
    };

    FullereneBatchView() = default;
    ~FullereneBatchView() = default;

    FullereneBatchView(const FullereneBatchView<T, K> &other) = default;
    FullereneBatchView(FullereneBatchView<T, K> &&other) = default;

    FullereneBatchView<T, K> &operator=(const FullereneBatchView<T, K> &other) = default;
    FullereneBatchView<T, K> &operator=(FullereneBatchView<T, K> &&other) = default;

    bool operator==(const FullereneBatchView<T, K> &other) const {d_ == other.d_ && m_ == other.m_;}

    FullereneBatchView(const FullereneBatchView<T, K> &other, size_t offset, int count) : 
        N_(other.N_), Nf_(other.Nf_),
        size_(count < 0 ? other.size() - offset : count){
        if ((offset + size_) > other.size()) { std::cerr << "offset is " << offset << " count is " << count << " size is " << other.size() << std::endl; throw std::out_of_range("Offset out of range");}
        fullerene_detail::construct_spans(offset, size_, d_.get_size_factors(other.N_, 1), d_.to_tuple(), other.d_.to_tuple());
        fullerene_detail::construct_spans(offset, size_, m_.get_size_factors(1), m_.to_tuple(), other.m_.to_tuple());
    }

    FullereneBatchView(const FullereneBatch<T, K> &batch, size_t offset = 0, int count = -1) : 
        N_(batch.N_), Nf_(batch.Nf_),
        size_(count < 0 ? batch.size() - offset : count){
        if ((offset + size_) > batch.capacity()) {std::cerr << "batch capacity is " << batch.capacity() << " offset is " << offset << " count is " << count << std::endl; throw std::out_of_range("Offset out of range");}
        fullerene_detail::construct_spans(offset, size_, d_.get_size_factors(batch.N_, 1), d_.to_tuple(), batch.d_.to_tuple());
        fullerene_detail::construct_spans(offset, size_, m_.get_size_factors(1), m_.to_tuple(), batch.m_.to_tuple());
    }

    FullereneBatchView(const Iterator& begin, const Iterator& end) : 
        N_(begin.get_view().N_), Nf_(begin.get_view().Nf_), size_(std::distance(begin, end)) {
        fullerene_detail::construct_spans(begin.get_index(), size_, d_.get_size_factors(N_, 1), d_.to_tuple(), begin.get_view().d_.to_tuple());
        fullerene_detail::construct_spans(begin.get_index(), size_, m_.get_size_factors(1), m_.to_tuple(), begin.get_view().m_.to_tuple());
    }

    Fullerene<T, K> operator[](size_t index) const {
        if (index >= size_) {assert(!"Index out of bounds");}
        FullereneDataMembers<Span, T, K> dnew;

        auto metatuple = m_.to_tuple();
        auto mnew = [&]<std::size_t... Is>(std::index_sequence<Is...>, auto& tup, std::size_t index) {
            return FullereneMetaMembers<ReferenceWrapper, K>{
                std::ref(std::get<Is>(tup)[index])...
            };
        }(std::make_index_sequence<std::tuple_size_v<decltype(metatuple)>>{}, metatuple, index);
        fullerene_detail::construct_spans(index, 1, d_.get_size_factors(N_, 1), dnew.to_tuple(), d_.to_tuple());
        return Fullerene<T, K>(dnew, mnew, N_, Nf_);
    }
                                
    // Getters for private members that are not supposed to be modified outside the class
    constexpr inline int size() const {return size_;}                            // Returns the size of the view
    constexpr inline int capacity() const {return size_;}                        // Returns the capacity of the view    
    constexpr inline Fullerene<T, K> front() const { return (*this)[0]; }         // Returns the first fullerene in the view
    constexpr inline Fullerene<T, K> back() const { return  (*this)[size() - 1]; } // Returns the last fullerene in the view
    
    constexpr inline auto begin() const { return Iterator(*this, 0); }
    constexpr inline auto end() const { return Iterator(*this, size_); }


    FullereneDataMembers<Span, T, K> d_;
    FullereneMetaMembers<Span, K> m_;

    const size_t N_;                      // Number of vertices in the cubic graph {1}
    const size_t Nf_;                     // Number of faces in the dual graph {1}
    const int size_; // Number of isomers in the view {1}

    //Output stream operator for FullereneBatchView
    template <typename U, typename V>
    friend std::ostream& operator<<(std::ostream& os, const FullereneBatchView<U, V>& batch);

    
};

namespace std {
    template <>
    struct iterator_traits<FullereneBatchView<>::Iterator> {
        using difference_type = std::ptrdiff_t;
        using value_type = typename FullereneBatchView<>::Iterator::value_type;
        using pointer = typename FullereneBatchView<>::Iterator::pointer;
        using reference = typename FullereneBatchView<>::Iterator::reference;
        using iterator_category = std::input_iterator_tag; // or the appropriate category
    };
}