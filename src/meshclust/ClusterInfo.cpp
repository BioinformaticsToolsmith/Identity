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
 * ClusterInfo.cpp
 *
 *  Created on: Jun 9, 2021
 *      Author: Hani Z. Girgis, PhD
 */
#include "ClusterInfo.h"

ClusterInfo::ClusterInfo(int id) {
	memberList = new std::vector<info>();
	repInfo = nullptr;
	repId = 0.0;
	repIndex = 0;
	identifier = id;
}

ClusterInfo::~ClusterInfo() {
	for (auto p : *memberList) {
		delete p.headerPtr;
	}
	delete memberList;
}

/*
 * member: the info of a sequence assigned to this cluster
 * idC: the identity score of the member with the mean of the cluster
 * idN: the identity score of the member with the mean of the closest cluster
 * membership: member (M), extended (E), outside (O)
 */
void ClusterInfo::addMember(std::string *memberInfo, double idC, double idN,
		int membership) {
	memberList->push_back( { memberInfo, idC, idN, membership });
	if (idC > repId) {
		repId = idC;
		repInfo = memberInfo;
		repIndex = memberList->size() - 1;
	}
}

void ClusterInfo::updateScoreWithCenter(double *list, int size) {
	if (size != memberList->size()) {
		std::cerr << "The size of the new list is incorrect." << std::endl;
		throw std::exception();
	}

	for (int i = 0; i < size; i++) {
		memberList->at(i).scoreWithCenter = list[i];
	}
}

void ClusterInfo::updateScoreWithNeighbor(double *list, int size) {
	if (size != memberList->size()) {
		std::cerr << "The size of the new list is incorrect." << std::endl;
		throw std::exception();
	}

	for (int i = 0; i < size; i++) {
		memberList->at(i).scoreWithNeighbor = list[i];
	}
}

std::vector<std::string*> ClusterInfo::getMemberList() {
	int size = memberList->size();
	std::vector<std::string*> v(size, nullptr);
	for (int i = 0; i < size; i++) {
		v[i] = memberList->at(i).headerPtr;
	}
	return v;
}

std::string* ClusterInfo::getCenter() {
	return memberList->at(repIndex).headerPtr;
}

int ClusterInfo::getRepIndex() {
	if (memberList->empty()) {
		std::cerr << "This cluster does not have any members and " << std::endl;
		std::cerr << "it does not have a representative sequence." << std::endl;
		throw std::exception();
	}
	return repIndex;
}

std::string ClusterInfo::toString() {
	std::ostringstream os;
	int l = memberList->size();
	for (int j = 0; j < l; j++) {
		auto p = memberList->at(j);
		char rounded[7];
		sprintf(rounded, "%.4f", p.scoreWithCenter);

		os << identifier << "\t" << *p.headerPtr << "\t"
				<< std::string(rounded);
		if (j == repIndex) {
			// This is the center
			os << "\tC";
		} else if (p.membership == ClusteringUtil::MEMBER) {
			// This sequence is within the threshold
			os << "\tM";
		} else if (p.membership == ClusteringUtil::EXTENDED) {
			// This sequence is within the relaxes threshold
			os << "\tE";
		} else if (p.membership == ClusteringUtil::OUTSIDE) {
			// This sequence is outside the relaxed threshold
			// However, the center of this cluster is the closest one to it
			os << "\tO";
		} else {
			std::cerr << "Invalid cluster membership: " << p.membership;
			std::cerr << std::endl;
			throw std::exception();
		}
		os << std::endl;
	}
	return os.str();
}

int ClusterInfo::getSize() {
	return memberList->size();
}

double ClusterInfo::silhoutte() {
	double s = 0.0;
	for (auto member : *memberList) {
		double distToC = 1.0 - member.scoreWithCenter;
		double distToN = 1.0 - member.scoreWithNeighbor;
		double m = distToN > distToC ? distToN : distToC;
		s += (distToN - distToC) / m;
	}
	return s;
}

double ClusterInfo::intra() {
	double s = 0.0;
	for (auto m : *memberList) {
		double distToC = 1.0 - m.scoreWithCenter;
		s += distToC;
	}

	if (memberList->empty()) {
		std::cerr << "An empty cluster cannot exist at this stage.";
		std::cerr << std::endl;
		throw std::exception();
	}

	return s / memberList->size();
}

void ClusterInfo::setIdentifier(int id) {
	identifier = id;
}

int ClusterInfo::getIdentifier() {
	return identifier;
}
