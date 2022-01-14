/*
 MeShClust v3.0 clusters sequences using the mean shift algorithm and alignment-free identity scores.

 Copyright (C) 2020-2022 Hani Z. Girgis, PhD

 Academic use: Affero General Public License version 1.

 Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 */

/*
 * ClusterEvaluator.h
 *
 *  Created on: Jun 10, 2021
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef SRC_MESHCLUST_CLUSTEREVALUATOR_H_
#define SRC_MESHCLUST_CLUSTEREVALUATOR_H_

#include <vector>
#include <iostream>
#include <math.h>

#include "ClusterInfo.h"
#include "../Matrix.h"

class ClusterEvaluator {
private:
	Matrix &ava;
	std::vector<ClusterInfo*> &clusterList;
	int total;

public:
	ClusterEvaluator(Matrix&, std::vector<ClusterInfo*>&, int);
	virtual ~ClusterEvaluator();
	double daviesBouldin();
	double dunn();
	double silhoutte();
	double clusterRatio();
	double intra();
	double inter();
	void printAll(bool isTab = false);
};

#endif /* SRC_MESHCLUST_CLUSTEREVALUATOR_H_ */
