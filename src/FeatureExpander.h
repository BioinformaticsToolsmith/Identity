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
 * FeatureExpander.h
 *
 *  Created on: Apr 14, 2020
 *      Author: Dr. Hani Z. Girgis
 */
#ifndef FEATUREEXPANDER_H_
#define FEATUREEXPANDER_H_

#include <vector>

#include "ITransformer.h"
#include "Matrix.h"
#include "Feature.h"
#include "FeatureSquared.h"
#include "FeaturePaired.h"

class FeatureExpander: public ITransformer {
private:
	/**
	 * This list includes pairs representing new squared and paired statistics
	 */
	int singleNum;

	std::vector<Feature*> fList;

public:
	FeatureExpander(const std::vector<Feature*>&);
	FeatureExpander(const Matrix&, const std::vector<Feature*>&);

	virtual ~FeatureExpander();

	/**
	 * Expands features in a matrix
	 */
	virtual Matrix transform(const Matrix&);

	virtual std::vector<Feature *> getFeatureList();
};

#endif /* FEATUREEXPANDER_H_ */
