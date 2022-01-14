/*
	Identity 2.0 calculates DNA sequence identity scores rapidly without alignment.

	Copyright (C) 2020-2022 Hani Z. Girgis, PhD

	Academic use: Affero General Public License version 1.

	Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

	Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
*/

/*
 * Node.cpp
 *
 *  Created on: Apr 28, 2020
 *      Author: Dr. Hani Zakaria Girgis
 *
 */

#include "Node.h"

Node::Node() {
	size = 0;
	list = nullptr;
}

/**
 * l: An array of features
 * s: Size of feature array
 */
Node::Node(int *l, int s) {
	list = l;
	size = s;

	// Pre-condition: l must be sorted
	for (int i = 0; i < size - 1; i++) {
		if (list[i] >= list[i + 1]) {
			std::cerr << "Node error: Unsorted array." << std::endl;
			std::cerr << *this << std::endl;
			throw std::exception();
		}
	}
}

Node::Node(const Node &n) {
	size = n.size;
	list = new int[size];
	for (int i = 0; i < size; i++) {
		list[i] = n.list[i];
	}
}

Node::Node(Node &&n) {
	size = n.size;
	list = n.list;
	n.list = nullptr;
}

Node::~Node() {
	if (list != nullptr) {
		delete[] list;
	}
}

Node& Node::operator=(const Node &n) {
	if (list != nullptr) {
		delete[] list;
	}

	size = n.size;
	list = new int[size];
	for (int i = 0; i < size; i++) {
		list[i] = n.list[i];
	}

	return *this;
}

/**
 * Construct a new node with the specified feature deleted
 * The resulting feature list is sorted
 */
Node Node::del(int f) {
	if (size == 0) {
		std::cerr << "Node error: Cannot perform delete on empty node.";
		std::cerr << std::endl;
		throw std::exception();
	}

	if (f < 0) {
		std::cerr << "Node error: A feature number must be non negative.";
		std::cerr << std::endl;
		throw std::exception();
	}

	int *l = new int[size - 1];
	int j = 0;
	for (int i = 0; i < size; i++) {
		if (list[i] != f) {
			l[j] = list[i];
			j++;
		}
	}

	return Node(l, size - 1);
}

/**
 * Construct a new node the specified feature added
 * The resulting feature list is sorted
 */
Node Node::add(int f) {
	if (f < 0) {
		std::cerr << "Node error: A feature number must be non negative.";
		std::cerr << std::endl;
		throw std::exception();
	}

	// Find index for new feature
	int i = 0;
	for (; i < size; i++) {
		if (list[i] > f) {
			break;
		}
	}
	// Construct the new node
	int *l = new int[size + 1];
	l[i] = f;

	for (int j = 0; j < i; j++) {
		l[j] = list[j];
	}

	for (int j = i + 1; j < size + 1; j++) {
		l[j] = list[j - 1];
	}

	return Node(l, size + 1);
}

/**
 * Perform the deletion and the addition operation on this node
 */
std::vector<Node> Node::expand(int fNum) {
	if (fNum <= 0) {
		std::cerr << "Node error: Feature number must be positive.";
		std::cerr << std::endl;
		throw std::exception();
	}

	std::vector<Node> l;
	l.reserve(fNum);
	for (int i = 0; i < size; i++) {
		l.push_back(this->del(list[i]));
	}

	int i = 0;
	int f = 0;
	while (f < fNum) {
		int limit = fNum;
		if (i < size) {
			limit = list[i];
			i++;
		}
		for (int j = f; j < limit; j++, f++) {
			Node n = this->add(j);
			l.push_back(n);
		}
		f++; // Skip list[i]
	}

	return l;
}

int Node::getSize() const {
	return size;
}

/**
 * Internal array must be sorted
 */
bool Node::operator==(const Node &other) const {
	bool r = false;
	if (other.getSize() == size) {
		const int *oList = other.getList();
		for (int i = 0; i < size; i++) {
			if (oList[i] != list[i]) {
				r = false;
				break;
			} else {
				r = true;
			}
		}
	}
	return r;
}

const int* Node::getList() const {
	return list;
}

/**
 * Credit: https://stackoverflow.com/questions/12840975/hashing-an-unordered-sequence-of-small-integers
 */
size_t NodeHasher::operator()(const Node &n) const {
	std::size_t result = 0;
	const int *l = n.getList();
	int s = n.getSize();

	for (int i = 0; i < s; i++) {
		result ^= l[i] + 0x9e3779b9 + (result << 6) + (result >> 2);
	}

	return result;
}

std::ostream& operator<<(std::ostream &os, const Node &n) {
	const int *l = n.getList();
	for (int i = 0; i < n.getSize(); i++) {
		os << l[i];
		if (i != n.getSize() - 1) {
			os << " ";
		}
	}
	return os;
}
