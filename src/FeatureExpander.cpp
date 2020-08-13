/*
	Identity calculates DNA sequence identity scores rapidly without alignment.

	Copyright (C) 2020 Hani Z. Girgis, PhD

	Academic use: Affero General Public License version 1.

	Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

	Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
*/

/*
 * FeatureExpander.cpp
 *
 *  Created on: Apr 14, 2020
 *      Author: Dr. Hani Z. Girgis
 *     Purpose: Copy original columns and add desired expanded ones (squared or paired).
 *
 */

#include "FeatureExpander.h"

/**
 * This constructor builds a list of all squared and paired statistics
 */
FeatureExpander::FeatureExpander(const Matrix &m,
		const std::vector<Feature*> &f) {

	singleNum = m.getNumCol();
	if (singleNum != f.size()) {
		std::cerr << "FeatureExpander error: " << std::endl;
		std::cerr << "Feature number does not match column number.";
		std::cerr << std::endl;
		throw std::exception();
	}

	int squaredNum = singleNum;
	int singleAndSquaredNum = singleNum + squaredNum;
	int pairedNum = (singleAndSquaredNum * (singleAndSquaredNum - 1)) / 2.0;

	fList = std::vector<Feature*>(singleNum + squaredNum + pairedNum, nullptr);
	// Copy single features
	int i = 0;
	for (i = 0; i < singleNum; i++) {
		fList[i] = new Feature(*f[i]);
	}

	// Add squared statistics to the feature list
	for (int c = 0; c < singleNum; c++, i++) {
		fList[i] = new FeatureSquared(fList[c]);
	}

	// Add paired statistics to the expansion list
	for (int c1 = 0; c1 < singleAndSquaredNum - 1; c1++) {
		for (int c2 = c1 + 1; c2 < singleAndSquaredNum; c2++, i++) {
			fList[i] = new FeaturePaired(fList[c1], fList[c2]);
		}
	}

	for (int h = 0; h < fList.size(); h++) {
		fList.at(h)->setTableIndex(h);
	}
}

/**
 * Constructor for optimized features
 */
FeatureExpander::FeatureExpander(const std::vector<Feature*> &f) {
	fList = copy(f);

	if (f.size() == 0) {
		std::cerr << "FeatureExpander error: No features are provided.";
		std::cerr << std::endl;
		throw std::exception();
	}

	singleNum = 0;
	for (auto ptr : fList) {
		if (ptr->getNumOfComp() == 0) {
			singleNum++;
		}
	}

	if (singleNum == 0) {
		std::cerr << "FeatureExpander error: No single features.";
		std::cerr << std::endl;
		throw std::exception();
	}

	// Test delete
	for (int i = singleNum; i < f.size(); i++) {
		auto feat = f.at(i);

		if (feat->getNumOfComp() == 1) {
			int c1 = feat->getCompOneIndex();
		}

		if (feat->getNumOfComp() == 2) {
			int c1 = feat->getCompOneIndex();
			int c2 = feat->getCompTwoIndex();
		}
	}
}

FeatureExpander::~FeatureExpander() {
	for (auto ptr : fList) {
		delete ptr;
	}
	fList.clear();
}

/**
 * Copy singles and expands features in a matrix
 */
Matrix FeatureExpander::transform(const Matrix &m) {
	if (m.getNumCol() != singleNum) {
		std::cerr << "FeatureExpander error: " << std::endl;
		std::cerr << "Column numbers do not match." << std::endl;
		throw std::exception();
	}

	int rowNum = m.getNumRow();
	int colNum = fList.size();
	Matrix t(rowNum, colNum);

	for (int i = 0; i < colNum; i++) {
		auto f = fList[i];

		if (i != f->getTableIndex()) {
			std::cerr << "FeatureExpander error: " << std::endl;
			std::cerr << "Index does not match feature table index.";
			std::cerr << std::endl;
			throw std::exception();
		}

		int n = f->getNumOfComp();
		// Handle singles
		if (n == 0) {
			for (int r = 0; r < rowNum; r++) {
				t(r, i) = m(r, i);
			}
		}
		// Handle squared
		else if (n == 1) {
			int c = f->getCompOneIndex();
			for (int r = 0; r < rowNum; r++) {
				double temp = m(r, c);
				t(r, i) = temp * temp;
			}
		}
		// Handle paired
		else if (n == 2) {
			int c1 = f->getCompOneIndex();
			int c2 = f->getCompTwoIndex();
			for (int r = 0; r < rowNum; r++) {
				t(r, i) = t(r, c1) * t(r, c2);
			}
		}
		// Handle unknown type
		else {
			std::cerr << "FeatureExpander error: " << std::endl;
			std::cerr << "Unexpected component number: " << n;
			std::cerr << std::endl;
			throw std::exception();
		}
	}

	return t;
}

std::vector<Feature*> FeatureExpander::getFeatureList() {
	return copy(fList);
}
