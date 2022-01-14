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
 * Cluster.h
 *
 *  Created on: Dec 23, 2020
 *      Author: Hani Z. Girgis, PhD
 *
 */

#ifndef SRC_MESHCLUST_CLUSTER_H_
#define SRC_MESHCLUST_CLUSTER_H_

#include <numeric>
#include <cstring>
#include <unordered_set>
#include <vector>

#include "../Util.h"

using namespace std;

template<class V>
class Cluster {
private:
	// Data info: All sequences in a set
	V **kHistList;
	uint64_t **monoHistList;

	int listSize;
	int kHistSize;
	int monoHistSize;

	// A list of identity scores between the mean and every sequence in a set
	double *identityList;

	// Those are the means that will be updated
	V *kHistMean;
	uint64_t *monoHistMean;

	// Those are the means of old members
	V *kHistOld;
	uint64_t *monoHistOld;

	// Number of points contributed to the mean of this cluster from
	// this block and all previous blocks
	// Updated in the merge step
	int contribution;

	// This is the number of members that contributed to the old members
	// Used in the merge and the update-accumulated-mean steps
	int oldContribution;

	// Number of points assigned to this cluster now or in the past
	int assignment;

	// Threshold score
	double threshold;

	// True if the identity list is up to date
	// This list needs updating when the representative sequence has changed
	bool isIdListUpToDate;

	// True if the mean has shifted
	bool hasShifted;

	// True if the id list is a row in the all-vs-all matrix
	bool isOriginalIdList;

	// This list is needed for selecting a representative sequence
	std::vector<int> * memberList;

public:
	Cluster(V**, uint64_t**, double*, int, int, int, double, int);
	Cluster(V**, uint64_t**, double*, int, int, int, double, V *kHist,
			uint64_t*, int, int);
	virtual ~Cluster();

	void init(V**, uint64_t**, double*, int, int, int, double);
	void setRepresentative(V*, uint64_t*, bool);
	void updateReferenceData(V **kList, uint64_t **monoList, int lSize);
	void shiftWeighted();

	void mergeSimple(std::vector<Cluster*>&);
	void updateAccumulatedMean();

	V* getKHistMean() const;
	uint64_t* getMonoHistMean() const;
	V* getKHistOld() const;
	uint64_t* getMonoHistOld() const;

	void setIdentityList(double *identityList);
	const double* getIdentityList();
	int getContribution() const;
	int getOldContribution() const;

	/**
	 * Return an estimated length based on the magnitude of the mono histogram
	 */
	int getLength() const;
	int getOldLength() const;

	void incrementAssignment();
	int getAssignment();

	bool getIsIdListUpToDate() const;
	bool getHasShifted() const;

	void copyIdentityList();

	bool getIsOriginalIdList() const;

	vector<int> * getMemberList() const;
	void clearMemberList();
};

#include "Cluster.cpp"

#endif /* SRC_MESHCLUST_CLUSTER_H_ */
