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
 * ITransform.h
 *
 *  Created on: Apr 14, 2020
 *      Author: Dr. Hani Z. Girgis
 */

#ifndef ITRANSFORMER_H_
#define ITRANSFORMER_H_

#include <vector>
#include <iostream>

#include "Matrix.h"
#include "Feature.h"
#include "FeatureSquared.h"
#include "FeaturePaired.h"

class ITransformer {

public:
//	ITransformer();
	virtual ~ITransformer() {
	}

	/**
	 *
	 */
	virtual Matrix transform(const Matrix&) = 0;

	/**
	 * Return a vector of pointers pointing to features info used by
	 * this transformer
	 */
	virtual std::vector<Feature*> getFeatureList() = 0;

	std::vector<Feature*> copy(const std::vector<Feature*>&);
};

#endif /* ITRANSFORMER_H_ */
