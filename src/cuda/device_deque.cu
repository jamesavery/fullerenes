template <typename T>
struct CuDeque
{
private:

    device_node_t front, back, q_size, capacity;
    T* array;

public:
    //This is a custom implementation of a deque class that is meant to be used on the device. 
//The deque is implemented using a circular buffer in a shared memory array. 
//The class is used to store the work queue of the thread blocks, and the class provides the
//necessary functions to allow the work queue to be used as a work stealing queue. 
//The class is implemented as a template class to allow the user to specify the type of data that the
//deque will hold.
__device__ CuDeque(T* memory, const device_node_t capacity): array(memory), front(0), back(0), q_size(0), capacity(capacity) {}
    
    /**
     * @brief  Returns the size of the queue. 
     * @param  None
     * @retval Size of the queue.
     */
    __device__ device_node_t size(){return q_size;}

    /**
     * @brief  Returns true if the queue is empty and false otherwise. 
     * @param  None
     * @retval True if the queue is empty and false otherwise.
     */
    __device__ bool empty(){
        return q_size == 0;
    }

    /**
     * @brief  Returns true if the queue is full and false otherwise. 
     * @param  None
     * @retval True if the queue is full and false otherwise.
     */
    __device__ bool full(){
        return q_size == capacity;
    }
    
    /**
     * @brief  This function is used to pop the first element of the queue.
     * @param  None
     * @retval First element of the queue
     */
    __device__ T pop_front(){
        if (!empty()){
            T return_val = array[front];
            front = (front + 1) % capacity ;
            q_size--;
            return return_val;
        }
        assert(false);
        return T(); //Compiler wants a return statement
    }

    /**
     * @brief Returns the last element of the queue and removes it from the queue.
     * @param None
     * @return The last element of the queue
     */
    __device__ T pop_back(){
        if (!empty())
        {
            T return_val = array[back];
            back = back > 0 ? back-1 : capacity-1;
            q_size--;
            return return_val;
        }
        assert(false);
        return T(); //Compiler wants a return statement
    }

    /** @brief Insert a value into the back of the queue
     *  @param val the value to insert
     */
    __device__ void push_back(T val){
        assert(!full());
        back = (back + 1) % capacity;
        array[back] = val;
        q_size++;
    }

    /** @brief Insert a value into the front of the queue
     *  @param val the value to insert
     */
    __device__ void push_front(T val){
        assert(!full());
        front = front > 0 ? front-1 : capacity-1;
        array[front] = val;
        q_size++;
    }
};