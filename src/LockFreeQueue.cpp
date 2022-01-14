/*
 * LockFreeQueue.cpp
 *
 *  Created on: Jul 7, 2020
 *      Author: Viktor Sehr and Bjorn Andrist (C++ High Performance)
 *      Edited by Hani Z. Girgis, PhD
 */

template<class T, size_t N>
LockFreeQueue<T,N>::LockFreeQueue() {
	assert(size_.is_lock_free());
}

template<class T, size_t N>
size_t LockFreeQueue<T,N>::size() const {
	return size_.load();
}

template<class T, size_t N>
void LockFreeQueue<T,N>::push(const T& t){
	if(size_.load() >= N){
		std::cerr << size_.load() << std::endl;
		throw std::overflow_error("LockFreeQueue error: Queue is full.");
	}
	buffer_[write_pos_] = t;
	write_pos_ = (write_pos_ + 1) % N;
	size_.fetch_add(1);
}

template<class T, size_t N>
T& LockFreeQueue<T,N>::front (){
	if(size_.load() == 0){
		throw std::underflow_error("LockFreeQueue error: Cannot call front on an empty queue.");
	}
	return buffer_[read_pos_];
}

template<class T, size_t N>
void LockFreeQueue<T,N>::pop (){
	if(size_.load() == 0){
		throw std::underflow_error("LockFreeQueue error: Cannot call pop on an empty queue.");
	}
	read_pos_ = (read_pos_ + 1) % N;
	size_.fetch_sub(1);
}
