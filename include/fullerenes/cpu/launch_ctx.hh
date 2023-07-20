#pragma once
#include "chrono"


  class LaunchCtx
  {
  public:
    LaunchCtx(){}
    LaunchCtx(int device) {}

    LaunchCtx& operator=(const LaunchCtx& other) {}
    ~LaunchCtx() {}
    int get_device_id() const {return 0;}
    //Queries the stream about its' work status.
    bool is_finished() const { return true; }
    //Returns true if the stream associated with the LaunchCtx is the default (NULL) stream and false otherwise.
    bool is_default_stream() const { return true; } 
    //Synchronizes this stream context with the calling host thread.
    void wait() const { }
    void start_timer() { }                	// TODO: Add timer
    std::chrono::nanoseconds stop_timer() {}	// TODO: Add timer

    //Synchronizes all streams, expensive operation, don't use unless necessary.
    static void wait_all() {}
    static int get_device_count() {return 1;}
    static int get_stream_count() {return 1;}
    static void clear_allocations() {}
  };


