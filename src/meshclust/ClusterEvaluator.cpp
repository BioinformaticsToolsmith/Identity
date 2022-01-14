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
 * ClusterEvaluator.cpp
 *
 *  Created on: Jun 10, 2021
 *      Author: Hani Z. Girgis, PhD
 */

#include "ClusterEvaluator.h"

ClusterEvaluator::ClusterEvaluator(Matrix &a, std::vector<ClusterInfo*> &l,
		int t) :
		ava(a), clusterList(l), total(t) {

}

ClusterEvaluator::~ClusterEvaluator() {
	// TODO Auto-generated destructor stub
}

double ClusterEvaluator::daviesBouldin() {
	int clustNum = clusterList.size();
	double sum = 0.0;

	for (int i = 0; i < clustNum; i++) {
		double max = -10000;
		auto clusterI = clusterList[i];
		for (int j = 0; j < clustNum; j++) {
			if (i != j) {
				auto clusterJ = clusterList[j];

				if (ava(i, j) >= 1.0) {
					std::cerr
							<< "Cannot have two identical clusters at this stage: "
							<< i << " " << j << " with score of " << ava(i, j)
							<< std::endl;
					throw std::exception();
				}

				double d = (clusterI->intra() + clusterJ->intra())
						/ (1.0 - ava(i, j));
				if (d > max) {
					max = d;
				}
			}
		}
		sum += max;
	}

	double result = std::numeric_limits<double>::infinity();
	if (clustNum > 0) {
		result = sum / clustNum;
	}

	return result;
}

double ClusterEvaluator::dunn() {
	// Find minimum inter
	double maxScore = -1;
	int rowNum = ava.getNumRow();
	int colNum = ava.getNumCol();
	for (int r = 0; r < rowNum; r++) {
		for (int c = 0; c < colNum; c++) {
			if (r != c && ava(r, c) > maxScore) {
				maxScore = ava(r, c);
			}
		}
	}
	double minInter = 1.0 - maxScore;

	// Find maximum intra
	double maxIntra = -1;
	for (auto clusterPtr : clusterList) {
		double intra = clusterPtr->intra();
		if (intra > maxIntra) {
			maxIntra = intra;
		}
	}
	double result = std::numeric_limits<double>::infinity();
	if (maxIntra > 0) {
		result = minInter / maxIntra;
	}

	return result;
}

double ClusterEvaluator::silhoutte() {
	double s = 0.0;
	int c = 0;
	for (auto clusterPtr : clusterList) {
		s += clusterPtr->silhoutte();
		c += clusterPtr->getSize();
	}

	double result = std::numeric_limits<double>::infinity();
	if (c > 0) {
		result = s / c;
	}

	return result;
}

/**
 * Return the ratio of clustered data to the total size
 */
double ClusterEvaluator::clusterRatio() {
	double clustered = 0.0;
	for (auto clusterPtr : clusterList) {
		clustered += clusterPtr->getSize();
	}
	return clustered / total;
}

double ClusterEvaluator::intra() {
	double sum = 0.0;
	for (auto clusterPtr : clusterList) {
		sum += clusterPtr->intra();
	}
	return 1.0 - (sum / clusterList.size());
}

/**
 * The inter-cluster identity of a cluster is the identity score to the closest one.
 * The returned score is the average over all clusters.
 */
double ClusterEvaluator::inter() {
	int rowNum = ava.getNumRow();
	int colNum = ava.getNumCol();
	double sum = 0.0;
	for (int r = 0; r < rowNum; r++) {
		double maxId = -1.0;
		for (int c = 0; c < colNum; c++) {
			if (r != c && ava(r, c) > maxId) {
				maxId = ava(r, c);
			}
		}

		if (maxId < 0.0 || maxId > 1.0) {
			std::cerr
					<< "Cannot determine the inter-cluster score. Max identity = "
					<< maxId << std::endl;
			ava.print();
			throw std::exception();
		}

		sum += maxId;
	}
	return sum / rowNum;
}

void ClusterEvaluator::printAll(bool isTab) {
	double db = daviesBouldin();
	double dn = dunn();
	double sil = silhoutte();
	double ratio = clusterRatio();
	double intraValue = intra();
	double interValue = inter();

	double quality = pow(
			(1 / db) * dn * ((1 + sil) / 2.0) * intraValue * (1 - interValue),
			1 / 5.0);

	std::cout << std::setprecision(4);

	if (isTab) {
		std::cout << "\t";
	}
	std::cout << "Davies-Bouldin index (lower is better): " << db << std::endl;
	if (isTab) {
		std::cout << "\t";
	}
	std::cout << "Dunn index (higher is better): " << dn << std::endl;
	if (isTab) {
		std::cout << "\t";
	}
	std::cout << "Silhoutte (higher is better): " << sil << std::endl;
	if (isTab) {
		std::cout << "\t";
	}
	std::cout << "Intra (higher is better): " << intraValue << std::endl;
	if (isTab) {
		std::cout << "\t";
	}
	std::cout << "Inter (lower is better): " << interValue << std::endl;
	if (isTab) {
		std::cout << "\t";
	}
	std::cout << "Cluster quality (higher is better): " << quality << std::endl;
	if (isTab) {
		std::cout << "\t";
	}
	std::cout << "Coverage (higher is better): " << ratio << std::endl;

}
