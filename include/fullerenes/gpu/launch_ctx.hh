#pragma once
#include "cuda_runtime_api.h"
#include "unordered_map"
#include "chrono"

class LaunchCtx
{
    public:
        LaunchCtx();
        LaunchCtx(int device);

        LaunchCtx& operator=(const LaunchCtx& other);
        ~LaunchCtx();
        int get_device_id() const;
        //Queries the stream about its' work status.
        bool is_finished() const;
        //Returns true if the stream associated with the LaunchCtx is the default (NULL) stream and false otherwise.
        bool is_default_stream() const;
        //Synchronizes this stream context with the calling host thread.
        void wait() const;
        void start_timer();
        std::chrono::nanoseconds stop_timer();
        //Synchronizes all streams, expensive operation, don't use unless necessary.
        static void wait_all();
        static int get_device_count();
        static int get_stream_count();
        static void clear_allocations();
        cudaStream_t stream;
    private:
        //The reason for pointer-pointer semantics is that the adress of the stream may change when the stream is created or destroyed.
        inline static std::unordered_map<int,cudaStream_t**> m_all_streams;
        //Increments every time a LaunchCtx is created, used as a way to uniquely identify streams, will break if more than 2147438647 streams are created.
        inline static int m_object_counter{};
        
        cudaEvent_t m_start;
        cudaEvent_t m_stop;
        int m_device_count;
        int m_device_id{};
        int m_unique_stream_idx{};
        inline static bool default_ctx_created{false};

};

