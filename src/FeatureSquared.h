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
 * FeatureSquared.h
 *
 *  Created on: May 20, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef FEATURESQUARED_H_
#define FEATURESQUARED_H_

#include <iostream>
#include "Feature.h"

class FeatureSquared: public Feature {
private:
	Feature *single;

public:
	FeatureSquared(Feature*);
	FeatureSquared();
	virtual ~FeatureSquared();
	FeatureSquared(const FeatureSquared &other);

	virtual const int getFunIndex() const;

	virtual int getNumOfComp();
	virtual int getCompOneIndex();
	virtual int getCompTwoIndex();

	virtual void setIsNeeded();
	virtual void setIsSelected();

	virtual void setCompOne(Feature*);
	virtual void setCompTwo(Feature*);
};

#endif /* FEATURESQUARED_H_ */
