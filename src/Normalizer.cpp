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
 * Normalizer.cpp
 *     Purpose: Normalize data between 0 and 1.
 *  Created on: Apr 4, 2020
 *      Author: Dr. Hani Z. Girgis
 */

#include "Normalizer.h"

/**
 * Normalizers (0--1) marked columns --- not all --- according to info in f
 */
Normalizer::Normalizer(const Matrix &m, const std::vector<Feature*> &f) :
		numCol(m.getNumCol()) {

	int numRow = m.getNumRow();
	// Pre conditions
	if (numRow == 0) {
		std::cerr << "Normalizer error: Matrix with no rows." << std::endl;
		throw std::exception();
	}

	if (numCol == 0) {
		std::cerr << "Normalizer error: Matrix with no columns." << std::endl;
		throw std::exception();
	}

	if (numCol != f.size()) {
		std::cerr << "Normalizer error: " << std::endl;
		std::cerr << "Feature number does not match column number.";
		std::cerr << std::endl;
		throw std::exception();
	}

	// Make a deep copy of feature list
	fList = copy(f);

	// Find the minimum and the maximum of each column
	for (int c = 0; c < numCol; c++) {
		auto f = fList[c];

		if (f->getIsNormalized()) {
			continue;
		}

		double min = m(0, c);
		double max = min;

		for (int r = 1; r < numRow; r++) {
			double i = m(r, c);
			if (i < min) {
				min = i;
			}
			if (i > max) {
				max = i;
			}
		}

		f->setNormP1(min);
		f->setNormP2(max);

		// A post condition
		if (max - min < 0.0) {
			std::cerr << "Normalizer error: Max - Min cannot be negative.";
			std::cerr << std::endl;
			throw std::exception();
		}
	}
}

Normalizer::Normalizer(const std::vector<Feature*> &f) :
		numCol(f.size()) {
	if (numCol == 0) {
		std::cerr << "Normalizer error: Matrix with no columns." << std::endl;
		throw std::exception();
	}

	// Make a deep copy of feature list
	fList = copy(f);
}

Normalizer::~Normalizer() {
	for (auto ptr : fList) {
		delete ptr;
	}
	fList.clear();
}

Matrix Normalizer::transform(const Matrix &m) {
	// Pre-condition
	if (m.getNumCol() != numCol) {
		std::cerr << "Normalizer error: Column numbers do not match.";
		std::cerr << "Expecting " << numCol << " but received "
				<< m.getNumCol();
		std::cerr << std::endl;
		throw std::exception();
	}

	int numRow = m.getNumRow();
	Matrix t(m);
	for (int c = 0; c < numCol; c++) {
		auto f = fList[c];
		if (!f->getIsNormalized()) {
			double min = f->getNormP1();
			double dist = f->getNormP2() - min;
			if (!Util::isEqual(dist, 0.0)) {
				for (int r = 0; r < numRow; r++) {
					double res = (m(r, c) - min) / dist;
					t(r, c) = res;
					// Trim if below 0 or above 1
					if (res < 0.0) {
						t(r, c) = 0.0;
					}
					if (res > 1.0) {
						t(r, c) = 1.0;
					}
				}
			} else {
				for (int r = 0; r < numRow; r++) {
					t(r, c) = 0;
				}
				std::cerr << "Normalizer warning: ";
				std::cerr << "Column " << c << " has the same value. ";
				std::cerr << "Values are already set to zero." << std::endl;
			}
		}
	}
	return t;
}

/**
 * Print minimums and maximums
 */
void Normalizer::printMinMaxLists() {
	std::cout << "Min\tMax" << std::endl;
	for (int i = 0; i < numCol; i++) {
		std::cout << i << "\t" << fList[i]->getNormP1() << "\t";
		std::cout << fList[i]->getNormP2();
		std::cout << std::endl;
	}
}

/**
 * Update features that have not been normalized yet
 * Update min (normP1), max (normP2), and isNormalized
 */
std::vector<Feature*> Normalizer::getFeatureList() {
	auto temp = copy(fList);
	for (auto f : temp) {
		f->setIsNormalized(true);
	}
	return temp;
}
