/*
 * LockFreeQueue.h
 *
 *  Created on: Jul 7, 2020
 *      Author: Viktor Sehr and Bjorn Andrist (C++ High Performance)
 *      Edited by Hani Z. Girgis, PhD
 */

#ifndef SRC_LOCKFREEQUEUE_H_
#define SRC_LOCKFREEQUEUE_H_

#include <cstdlib>
#include <atomic>
#include <array>
#include <assert.h>

template<class T, size_t N>
class LockFreeQueue {
private:
	int read_pos_ = 0;
	int write_pos_ = 0;
	std::atomic<size_t> size_ { 0 };
	std::array<T, N> buffer_ { };

public:
	LockFreeQueue();
	size_t size() const;
	void push(const T &t);
	T& front();
	void pop();
};

#include "LockFreeQueue.cpp"

#endif /* SRC_LOCKFREEQUEUE_H_ */
