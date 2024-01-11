template <typename T>
struct CuDeque
{
private:

    int front, back, q_size, capacity;
    T* array;

public:
    //This is a custom implementation of a deque class that is meant to be used on the device. 
//The deque is implemented using a circular buffer in a shared memory array. 
//The class is used to store the work queue of the thread blocks, and the class provides the
//necessary functions to allow the work queue to be used as a work stealing queue. 
//The class is implemented as a template class to allow the user to specify the type of data that the
//deque will hold.
__device__ CuDeque(T* memory, const int capacity): array(memory), front(-1), back(0), q_size(0), capacity(capacity) {}
    
    /**
     * @brief  Returns the size of the queue. 
     * @param  None
     * @retval Size of the queue.
     */
    __device__ int size(){return q_size;}

    /**
     * @brief  Returns true if the queue is empty and false otherwise. 
     * @param  None
     * @retval True if the queue is empty and false otherwise.
     */
    __device__ bool empty(){
        return (front == -1);
    }

    /**
     * @brief  Returns true if the queue is full and false otherwise. 
     * @param  None
     * @retval True if the queue is full and false otherwise.
     */
    __device__ bool full(){
        return (front == 0 && back == capacity-1) || (front == back+1);
    }
    
    /**
     * @brief  This function is used to pop the first element of the queue.
     * @param  None
     * @retval First element of the queue
     */
    __device__ T pop_front(){
        if (empty()){ assert(false); return T();} 
        T return_val = array[front];
        if(front == back) {
            front = -1;
            back = -1;
        } 
        else if (front == capacity-1) front = 0;
        else front = front+1;
        q_size--;
        return return_val;
    }

    /**
     * @brief Returns the last element of the queue and removes it from the queue.
     * @param None
     * @return The last element of the queue
     */
    __device__ T pop_back(){
        if (empty()){ assert(false); return T();}
        T return_val = array[back];
        if(front == back) {
            front = -1;
            back = -1;
        } 
        else if (back == 0) back = capacity-1;
        else back = back-1;
        q_size--;
        return return_val;
    }

    /** @brief Insert a value into the back of the queue
     *  @param val the value to insert
     */
    __device__ void push_back(T val){
        assert(!full());
        if (front == -1) {
            front = 0;
            back = 0;
        }
        else if (back == capacity-1) back = 0;
        else back = back+1;
        array[back] = val;
        q_size++;
    }

    /** @brief Insert a value into the front of the queue
     *  @param val the value to insert
     */
    __device__ void push_front(T val){
        assert(!full());
        if (front == -1) {
            front = 0;
            back = 0;
        }
        else if (front == 0) front = capacity-1;
        else front = front-1;
        array[front] = val;
        q_size++;
    }
};