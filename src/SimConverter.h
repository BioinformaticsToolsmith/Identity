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
 * SimConverter.h
 *
 *  Created on: Apr 14, 2020
 *      Author: Dr. Hani Z. Girgis
 */
#ifndef SIMCONVERTER_H_
#define SIMCONVERTER_H_

#include "ITransformer.h"
#include "Matrix.h"
#include "Feature.h"

class SimConverter: public ITransformer {
private:
	const int numCol;
	// A vector of statistics numbers, i.e. statistic names
	std::vector<Feature*> fList;

public:
	SimConverter(const Matrix&, const std::vector<Feature*>&);
	SimConverter(const std::vector<Feature*>&);
	virtual ~SimConverter();

	/**
	 * Convert distances in matrix to similarities
	 */
	virtual Matrix transform(const Matrix&);

	/**
	 * Does not change it. It reads only.
	 * Throws exception if it is called.
	 */
	std::vector<Feature*> getFeatureList();
};

#endif /* SIMCONVERTER_H_ */
