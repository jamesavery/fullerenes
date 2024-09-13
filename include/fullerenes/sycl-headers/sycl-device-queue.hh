#pragma once
#include <memory>
#include <vector>


enum class Policy
{
    SYNC,
    ASYNC
};

enum class DeviceType
{
    CPU,
    GPU,
    ACCELERATOR,
    HOST,
    NUM_DEV_TYPES
};

enum class DeviceProperty
{
    MAX_WORK_GROUP_SIZE,
    MAX_CLOCK_FREQUENCY,
    MAX_COMPUTE_UNITS,
    MAX_MEM_ALLOC_SIZE,
    GLOBAL_MEM_SIZE, 
    LOCAL_MEM_SIZE,
    MAX_CONSTANT_ARGS,
    MAX_NUM_SUB_GROUPS,
    MAX_SUB_GROUP_SIZE,
    NUMBER_OF_PROPERTIES
};

struct Device{
    static std::vector<Device> get_devices(DeviceType type);

    Device() = default;

    Device(size_t idx, DeviceType type) : idx(idx), type(type) {}

    inline static std::vector<Device> cpus_ =         get_devices(DeviceType::CPU);
    inline static std::vector<Device> gpus_ =         get_devices(DeviceType::GPU);
    inline static std::vector<Device> accelerators_ = get_devices(DeviceType::ACCELERATOR);

    std::string get_name() const;
    size_t get_property(DeviceProperty property) const;
    

    size_t     idx  = 0;
    DeviceType type = DeviceType::HOST;
};

struct SyclEventImpl;

struct SyclEvent {
    std::unique_ptr<SyclEventImpl> impl_;

    SyclEvent();
    ~SyclEvent();
    SyclEvent& operator=(SyclEventImpl&& impl);
    SyclEvent::SyclEvent(SyclEventImpl&& impl);
    SyclEvent(SyclEvent&& other);
    SyclEvent& operator=(SyclEvent&& other);
    void wait() const;
    SyclEventImpl* operator->() const;
    SyclEventImpl& operator*() const;
};

struct SyclQueueImpl;

struct SyclQueue{

    /*  
        Default constructor and destructor must be declared ONLY here and defined in the implementation file. 
        This is necessary because SyclQueueImpl is an incomplete type in this header file.
    */
    SyclQueue(); 
    ~SyclQueue();

    SyclQueue(Device device, bool in_order = true);
    SyclQueue(SyclQueue&& other) = default;
    SyclQueue& operator=(SyclQueue&& other) = default;
    SyclQueue(const SyclQueue& other) = delete;
    SyclQueue& operator=(const SyclQueue& other) = delete;
    // Copy constructor and assignment operator are deleted because the unique_ptr is non-copyable


    SyclQueueImpl* operator->() const;
    SyclQueueImpl& operator*() const;
    
    void enqueue(SyclEvent& event);
    SyclEvent get_event();

    template <typename EventContainer>
    void enqueue(EventContainer& events){
        for(auto& event : events){
            enqueue(event);
        }
    }

    void wait() const;
    void wait_and_throw() const;
    
    std::unique_ptr<SyclQueueImpl> impl_;
    
    const Device device_;
    const bool in_order_;
};

