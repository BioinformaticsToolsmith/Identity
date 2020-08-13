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
 * SimConverter.cpp
 *
 *  Created on: Apr 14, 2020
 *      Author: Dr. Hani Z. Girgis
 *
 */

#include "SimConverter.h"

/**
 * Convert a normalized distance (0--1) to a similarity measure by
 * subtracting it from 1.
 * l: A list of features, each of which contains information
 */
SimConverter::SimConverter(const Matrix &m, const std::vector<Feature*> &l) :
		numCol(m.getNumCol()) {
	if (l.size() != numCol) {
		std::cerr << "SimConverter error: " << std::endl;
		std::cerr << "Number of features does not match number of columns.";
		std::cerr << std::endl;
		throw std::exception();
	}

	if (numCol == 0) {
		std::cerr << "SimConverter error: Matrix with no columns." << std::endl;
		throw std::exception();
	}

	fList = copy(l);
}

/**
 * Constructor for optimized features
 */
SimConverter::SimConverter(const std::vector<Feature*> &l) :
		numCol(l.size()) {

	if (numCol == 0) {
		std::cerr << "SimConverter error: Matrix with no columns." << std::endl;
		throw std::exception();
	}

	fList = copy(l);
}

SimConverter::~SimConverter() {
	for (auto ptr : fList) {
		delete ptr;
	}
	fList.clear();
}

/**
 * Convert distances in matrix to similarities
 */
Matrix SimConverter::transform(const Matrix &m) {
	if (m.getNumCol() != numCol) {
		std::cerr << "SimConverter error: " << std::endl;
		std::cerr << "Number of features does not match number of columns.";
		std::cerr << std::endl;
		throw std::exception();
	}

	int numRow = m.getNumRow();
	Matrix t(m);

	for (int col = 0; col < numCol; col++) {
		auto f = fList[col];
		if (f->getIsDistance() && !f->getIsConverted()) {
			for (int row = 0; row < numRow; row++) {
				t(row, col) = 1 - m(row, col);
			}
		}
	}

	return t;
}

std::vector<Feature*> SimConverter::getFeatureList() {
	auto temp = copy(fList);
	for (auto f : temp) {
		if (f->getIsDistance() && !f->getIsConverted()) {
			f->setIsConverted(true);
		}
	}
	return temp;
}
