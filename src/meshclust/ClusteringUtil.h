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
 * ClusteringUtil.h
 *
 *  Created on: March 6, 2021
 *      Author: Dr. Hani Z. Girgis
 */

#ifndef CLUSTERINGUTIL_H_
#define CLUSTERINGUTIL_H_

#include <vector>
#include <queue>
#include "../Matrix.h"

class ClusteringUtil {
public:

	static const int MEMBER = 0;
	static const int EXTENDED = 1;
	static const int OUTSIDE = 2;

	static std::pair<std::vector<int>, int> findConnectedComponents(Matrix &m,
			double threshold) {
		// Find connected components using breadth first search
		// Credit: Algorithm idea and pseudo code from stackoverflow.com
		int r = m.getNumRow();
		int c = m.getNumCol();
		if (r != c) {
			std::cerr << "Util::findConnectedComponents needs a squared "
					<< "matrix. Received " << r << " rows and " << c
					<< " columns." << std::endl;
			throw std::exception();
		}

		int compNum = 0;
		std::vector<int> flagList(r, 0);
		for (int i = 0; i < r; i++) {
			if (flagList[i] == 0) {
				compNum++;
				std::queue<int> q;
				q.push(i);
				bool inQueueList[r] { false };
				inQueueList[i] = true;
				while (!q.empty()) {
					int v = q.front();
					q.pop();
					flagList[v] = compNum;
					for (int j = v + 1; j < r; j++) {
						// ToDo: Change at to ()
						if (!inQueueList[j] && flagList[j] == 0
								&& m.at(v, j) >= threshold) {
							q.push(j);
							inQueueList[j] = true;
						}
					}
				}
			}
		}
		return std::make_pair(flagList, compNum);
	}
};

#endif /* CLUSTERINGUTIL_H_ */
