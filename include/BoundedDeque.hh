#include <vector>
#include <cstdint>
#include <iostream>

template <typename T>
class BoundedDeque{

    // A bounded deque implemented with a circular buffer
    // There are no safety checks for example for pushing more elements than
    // the maximum capacity, or accessing an empty deque, so be careful. 

    private:

    int64_t mod_increment(int64_t i){
        return (i + 1) % buf.size();
    }

    int64_t mod_decrement(int64_t i){
        return (i - 1 + buf.size()) % buf.size();
    }

    std::vector<T> buf; // Circular buffer
    int64_t front_idx; // Index where the next element added to the front of the queue would go to 
    int64_t back_idx; // Index where the next element added to the back to the queue would go to
    int64_t n_elements;

    // buf: -----XXXXX---
    //          ^     ^
    //        front  back

    public:

    BoundedDeque(int64_t max_size) : buf(max_size), front_idx(max_size-1), back_idx(0), n_elements(0) {}

    const T& back(){
        return buf[mod_decrement(back_idx)];
    }

    const T& front(){
        return buf[mod_increment(front_idx)];
    }

    int64_t size(){
        return n_elements;
    }

    void push_back(const T& x){
        buf[back_idx] = x;
        back_idx = mod_increment(back_idx);
        n_elements++;
    }

    void push_front(const T& x){
        buf[front_idx] = x;
        front_idx = mod_decrement(front_idx);
        n_elements++;
    }

    void pop_front(){
        front_idx = mod_increment(front_idx);
        n_elements--;
    }

    void pop_back(){
        back_idx = mod_decrement(back_idx);
        n_elements--;
    }

    void clear(){
        n_elements = 0;
        front_idx = buf.size()-1;
        back_idx = 0;
    }

};
