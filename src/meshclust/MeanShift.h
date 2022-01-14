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
 * MeanShift.h
 *
 *  Created on: Dec 23, 2020
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef SRC_MESHCLUST_MEANSHIFT_H_
#define SRC_MESHCLUST_MEANSHIFT_H_

#include <string>
#include <vector>

#include "Cluster.h"
#include "ClusteringUtil.h"
#include "../Parameters.h"
#include "../IdentityCalculator.h"
#include "../Matrix.h"

typedef std::vector<std::pair<std::string*, std::string*> > Block;

template<class V>
class MeanShift {
private:
	IdentityCalculator<V> &identity;
	int threadNum;
	// Applied at the shift and assign steps
	double threshold;
	// Applied only at the merge step
	double mergeThreshold;

	int size;
	V **kHistList;
	uint64_t **monoHistList;
	std::string **infoList;
	int *lenList;

	int *assignList = nullptr;

	bool isDataCleared = false;

	//double mergeThreshold;

	// A matrix of all-versus-all identity scores
	Matrix a;
	bool isMatrixCleared;
	vector<Cluster<V>*> *clusterList;
	void updateIdentityList();
	void initData(Block *block);
	void initClusters();

	void start();
	// Cluster<V>* mergeGreedyHelper(std::vector<Cluster<V>*>&);

	bool isLowIdentity;

public:
	MeanShift(Block*, IdentityCalculator<V>&, int, double threshold);

	MeanShift(std::tuple<V**, uint64_t**, std::string**, int*, int>,
			IdentityCalculator<V>&, int, double threshold);

	virtual ~MeanShift();
	void clearData();
	void run(int itrNum, bool canAssign);
	//void runOnce();
	void removeSingles();
	void removeEmpty();
	void updateReferenceData(Block*);
	void addClusters(const vector<Cluster<V>*> *oldClusterList);
	std::tuple<V**, uint64_t**, std::string**, int*, int> findUnassignedData();
	void shift();
	//void merge();
	void mergeGreedy();
	void selectRepresentative();
	void updateAccumulatedMean();
	void assign();
	void copyIdentityList();

	Matrix calcAllCenterVsAllCenter();

	const vector<Cluster<V>*>* getClusterList() const;
	//double getThreshold();
	//void checkClusters();
};

#include "MeanShift.cpp"

#endif /* SRC_MESHCLUST_MEANSHIFT_H_ */
