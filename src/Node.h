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
 * Node.h
 *
 *  Created on: Apr 28, 2020
 *      Author: Dr. Hani Z. Girgis
 *     Purpose: A graph node for the best-first search algorithm
 */

#ifndef NODE_H_
#define NODE_H_

#include <vector>
#include <iostream>
#include <algorithm>

class Node {
private:
	int size;
	int *list;

public:
	Node();
	Node(int*, int);
	Node(const Node&);
	Node(Node&&);
	virtual ~Node();

	Node& operator=(const Node &n);

	int getSize() const;
	bool operator==(const Node &other) const;
	const int* getList() const;

	Node del(int);
	Node add(int);
	std::vector<Node> expand(int);
};

class NodeHasher {
public:
	/**
	 * A hashing function for Node
	 */
	size_t operator()(const Node &p) const;
};

std::ostream& operator<<(std::ostream &os, const Node &n);

#endif /* NODE_H_ */
