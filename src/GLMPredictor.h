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
 * GLMPredictor.h
 *
 *  Created on: Jul 9, 2020
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef SRC_GLMPREDICTOR_H_
#define SRC_GLMPREDICTOR_H_

#include <vector>
#include <iostream>

#include "Feature.h"

class GLMPredictor {
private:
	bool canDelete; // True is arrays have been allocated on the heap

	bool isClassification;
	int singleFeatNum;
	int featNum;
	double bias;
	int distNum;
	int selectNum;

	double *minList;
	double *maxMinList; // max - min
	int *distIndexList;
	int *selectedIndexList;
	double *wList;
	std::pair<int, int> *expList;

	void copy(const GLMPredictor&);

public:
	GLMPredictor();
	GLMPredictor(std::vector<Feature*>&, bool);
	GLMPredictor(const GLMPredictor&);
	virtual ~GLMPredictor();

	GLMPredictor& operator=(const GLMPredictor&);

	inline double calculateIdentity(double *data) {
		// Normalize and trim singles
		for (int i = 0; i < singleFeatNum; i++) {
			double d = (data[i] - minList[i]) / maxMinList[i];
			if (d > 1.0) {
				d = 1.0;
			}
			if (d < 0.0) {
				d = 0.0;
			}
			data[i] = d;
		}

		// Convert distances to similarities
		for (int i = 0; i < distNum; i++) {
			int index = distIndexList[i];
			data[index] = 1 - data[index];
		}

		// Expand
		for (int i = singleFeatNum; i < featNum; i++) {
			auto p = expList[i];
			data[i] = data[p.first] * data[p.second];
		}

		// Normalize and trim squares and pairs
		for (int i = singleFeatNum; i < featNum; i++) {
			double d = (data[i] - minList[i]) / maxMinList[i];
			if (d > 1.0) {
				d = 1.0;
			}
			if (d < 0.0) {
				d = 0.0;
			}
			data[i] = d;
		}

		// Calculate identity
		double res = bias;
		for (int i = 0; i < selectNum; i++) {
			res += (wList[i] * data[selectedIndexList[i]]);
		}

		if (isClassification) {
			res = res >= 0.5 ? 1.0 : 0.0;
		}

		return res;
	}

	int getFeatNum() const;
	int* getFunIndexArray() const;
	int getSingleFeatNum() const;
};

#endif /* SRC_GLMPREDICTOR_H_ */
