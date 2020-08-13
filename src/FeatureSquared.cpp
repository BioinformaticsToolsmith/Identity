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
 * FeatureSquared.cpp
 *
 *  Created on: May 20, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#include "FeatureSquared.h"

FeatureSquared::FeatureSquared(Feature *s) :
		Feature(-1, s->getName() + "^2", false) {
	single = s;
}

FeatureSquared::FeatureSquared(const FeatureSquared &o) :
		Feature(o) {
	single = o.single;
}

FeatureSquared::~FeatureSquared() {}

const int FeatureSquared::getFunIndex() const {
	std::cerr << "FeatureSquared error: " << std::endl;
	std::cerr << "A squared feature does not have funIndex.";
	std::cerr << std::endl;
	throw std::exception();
}

int FeatureSquared::getNumOfComp() {
	return 1;
}

int FeatureSquared::getCompOneIndex() {
	return single->getTableIndex();
}

int FeatureSquared::getCompTwoIndex() {
	std::cerr << "FeatureSquared error: " << std::endl;
	std::cerr << "A squared feature does not have a second component.";
	std::cerr << std::endl;
	throw std::exception();
}

void FeatureSquared::setCompOne(Feature *f) {
	single = f;
}

void FeatureSquared::setCompTwo(Feature *f) {
	std::cerr << "FeatureSquared error: " << std::endl;
	std::cerr << "A squared feature does not have a second component.";
	std::cerr << std::endl;
	throw std::exception();
}

void FeatureSquared::setIsSelected() {
	isSelected = true;
	single->setIsNeeded();
}

void FeatureSquared::setIsNeeded() {
	isNeeded = true;
	single->setIsNeeded();
}
