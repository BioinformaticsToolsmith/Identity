/*
 Identity calculates DNA sequence identity scores rapidly without alignment.

 Copyright (C) 2020-2022 Hani Z. Girgis, PhD

 Academic use: Affero General Public License version 1.

 Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 */

/*
 * Reservoir.h
 *
 *  Created on: Jul 5, 2021
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef SRC_MESHCLUST_RESERVOIR_H_
#define SRC_MESHCLUST_RESERVOIR_H_

#include <string>
#include <vector>
#include <algorithm>

template<class V>
class Reservoir {
private:
	std::vector<V*> kHistList;
	std::vector<uint64_t*> monoHistList;
	std::vector<std::string*> infoList;
	std::vector<int> lenList;
	int seed;

public:
	Reservoir();
	virtual ~Reservoir();
	void add(std::tuple<V**, uint64_t**, std::string**, int*, int>);
	std::tuple<V**, uint64_t**, std::string**, int*, int> remove(int);
	int size();
	void shuffle();
};

#include "Reservoir.cpp"

#endif /* SRC_MESHCLUST_RESERVOIR_H_ */
