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
 * Evaluator.cpp
 *
 *  Created on: May 16, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#include "Evaluator.h"

/**
 * Calculates the accuracy
 * oLabels: Original labels
 * pLabels: Predicted labels
 */
double Evaluator::acc(const Matrix &oLabels, const Matrix &pLabels) {
	int r = oLabels.getNumRow();
	if (r != pLabels.getNumRow()) {
		std::cerr << "Evaluator error: Original and predicted ";
		std::cerr << "labels have different sizes." << std::endl;
		throw std::exception();
	}

	const double *o = oLabels.getArray();
	const double *p = pLabels.getArray();

	double s = 0;
	for (int i = 0; i < r; i++) {
		if (Util::isEqual(o[i], p[i])) {
			s++;
		}
	}
	return s / r;
}

/**
 * oLabels: Original labels
 * pLabels: Predicted labels
 * l: Label of interest
 */
double Evaluator::helper(const Matrix &oLabels, const Matrix &pLabels,
		double l) {
	int r = oLabels.getNumRow();
	if (r != pLabels.getNumRow()) {
		std::cerr << "Evaluator error: Original and predicted ";
		std::cerr << "labels have different sizes." << std::endl;
		throw std::exception();
	}

	const double *o = oLabels.getArray();
	const double *p = pLabels.getArray();

	double t = 0;
	double s = 0;
	for (int i = 0; i < r; i++) {
		if (Util::isEqual(o[i], l)) {
			t++;
			if (Util::isEqual(o[i], p[i])) {
				s++;
			}
		}
	}
	return s / t;
}

/**
 * Calculates the sensitivity
 * oLabels: Original labels
 * pLabels: Predicted labels
 */
double Evaluator::sens(const Matrix &oLabels, const Matrix &pLabels) {
	return helper(oLabels, pLabels, 1.0);
}

/**
 * Calculates the specificity
 * oLabels: Original labels
 * pLabels: Predicted labels
 */
double Evaluator::spec(const Matrix &oLabels, const Matrix &pLabels) {
	return helper(oLabels, pLabels, 0.0);
}

/**
 * Mean Absolute Error
 */
double Evaluator::mae(const Matrix &oLabels, const Matrix &pLabels) {
	int r = oLabels.getNumRow();
	if (r != pLabels.getNumRow()) {
		std::cerr << "Evaluator error: Original and predicted ";
		std::cerr << "labels have different sizes." << std::endl;
		throw std::exception();
	}

	const double *o = oLabels.getArray();
	const double *p = pLabels.getArray();

	double s = 0.0;
	for (int i = 0; i < r; i++) {
		s += fabs(o[i] - p[i]);
	}
	return s / r;
}

/**
 * Mean Squared Error
 */
double Evaluator::mse(const Matrix &oLabels, const Matrix &pLabels) {
	int r = oLabels.getNumRow();
	if (r != pLabels.getNumRow()) {
		std::cerr << "Evaluator error: Original and predicted ";
		std::cerr << "labels have different sizes." << std::endl;
		throw std::exception();
	}

	const double *o = oLabels.getArray();
	const double *p = pLabels.getArray();

	double s = 0.0;
	for (int i = 0; i < r; i++) {
		double temp = fabs(o[i] - p[i]);
		temp *= temp;
		s += temp;
	}
	return s / r;
}
