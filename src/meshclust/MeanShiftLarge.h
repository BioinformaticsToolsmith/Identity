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
 * MeanShiftLarge.h
 *
 *  Created on: Dec 23, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef SRC_MESHCLUST_MEANSHIFTLARGE_H_
#define SRC_MESHCLUST_MEANSHIFTLARGE_H_

#include <vector>
#include <fstream>

#include "MeanShift.h"
#include "ClusterInfo.h"
#include "ClusterEvaluator.h"
#include "Reservoir.h"
#include "ClusteringUtil.h"
#include "../Matrix.h"
#include "../FastaReader.h"
#include "../IdentityCalculator.h"

template<class V>
class MeanShiftLarge {
private:
	string dbFile;
	string outFile;
	int blockSize;
	int vBlockSize;
	int passNum;
	IdentityCalculator<V> &identity;
	int threadNum;
	double threshold;
	MeanShift<V> *ms;
	Reservoir<V> reservoir;
	bool canAssignAll;
	bool canEvaluate;
	bool canRelax;

	// Sequences processed so far
	long seqNum;

	void clusterReservoir();
	void assign();
	Matrix calcAllCenterVsAllCenter();

	int printCounter = 0;
	int assignCounter = 1;
	void printStatus(bool printAnyWay=false);

public:
	MeanShiftLarge(string, int, int, int, IdentityCalculator<V>&, int, bool, string,
			double t, bool, bool);
	virtual ~MeanShiftLarge();
};

#include "MeanShiftLarge.cpp"

#endif /* SRC_MESHCLUST_MEANSHIFTLARGE_H_ */
