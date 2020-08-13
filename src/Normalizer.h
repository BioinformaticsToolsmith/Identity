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
 * Normalizer.h
 *
 *  Created on: Apr 4, 2020
 *      Author: Dr. Hani Zakaria Girgis
 */

#ifndef NORMALIZER_H_
#define NORMALIZER_H_

#include <iostream>
#include <vector>

#include "ITransformer.h"
#include "Matrix.h"
#include "Util.h"
#include "Feature.h"

class Normalizer: public ITransformer {
private:
	// A list of features
	std::vector<Feature*> fList;

	// Number of required columns
	const int numCol;

public:
	Normalizer(const std::vector<Feature*>&);
	Normalizer(const Matrix&, const std::vector<Feature*>&);

	virtual ~Normalizer();

	/**
	 * Normalize items in matrix between 0 and 1.
	 */
	virtual Matrix transform(const Matrix&);

	/**
	 * Normalize items in vector between 0 and 1.
	 */
	// virtual std::vector<double> transform(const std::vector<double>&);

	/**
	 * Update min (normP1), max (normP2)
	 */
	virtual std::vector<Feature*> getFeatureList();

	/**
	 * Print minList and maxList
	 */
	void printMinMaxLists();
};

#endif /* NORMALIZER_H_ */
