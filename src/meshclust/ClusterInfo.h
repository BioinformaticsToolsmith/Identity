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
 * ClusterInfo.h
 *
 *  Created on: Jun 9, 2021
 *      Author: Hani Z. Girgis, PhD
 */
#ifndef SRC_MESHCLUST_CLUSTERINFO_H_
#define SRC_MESHCLUST_CLUSTERINFO_H_

#include <vector>
#include <string>
#include <iostream>

#include "ClusteringUtil.h"

typedef struct {
	std::string *headerPtr;
	double scoreWithCenter;
	double scoreWithNeighbor;
	int membership;
} info;

class ClusterInfo {
private:
	// Member list
	std::vector<info> *memberList;
	std::string *repInfo;
	double repId; // Identity score of the representative center with the mean
	int repIndex;
	int identifier;

public:
	ClusterInfo(int id = 0);
	virtual ~ClusterInfo();
	void addMember(std::string*, double, double idN, int membership);
	void updateScoreWithCenter(double*, int);
	void updateScoreWithNeighbor(double*, int);
	int getRepIndex();
	std::string toString();

	int getSize();
	double silhoutte();
	double intra();
	double daviesBouldin();
	std::vector<std::string*> getMemberList();
	std::string* getCenter();
	void setIdentifier(int);
	int getIdentifier();
};

#endif /* SRC_MESHCLUST_CLUSTERINFO_H_ */
