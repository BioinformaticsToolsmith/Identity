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
 * IPredictor.h
 *
 *  Created on: Mar 31, 2020
 *      Author: Dr. Hani Z. Girgis
 */
#include <tuple>
#include <math.h>

#include "Matrix.h"
#include "Util.h"

#ifndef EVALUATOR_H_
#define EVALUATOR_H_

class Evaluator {
public:
	static double acc(const Matrix &oLabels, const Matrix &pLabels);
	static double sens(const Matrix &oLabels, const Matrix &pLabels);
	static double spec(const Matrix &oLabels, const Matrix &pLabels);
	/**
	 * Helper function to sens and spec
	 */
	static double helper(const Matrix &oLabels, const Matrix &pLabels,
			double l);
	/**
	 * Mean Absolute Error
	 */
	static double mae(const Matrix &oLabels, const Matrix &pLabels);
	/**
	 * Mean Squared Error
	 */
	static double mse(const Matrix &oLabels, const Matrix &pLabels);

public:
	virtual ~Evaluator() {}
};

#endif /* EVALUATOR_H_ */
