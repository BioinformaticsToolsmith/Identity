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
 * FeaturePaired.cpp
 *
 *  Created on: May 20, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#include "FeaturePaired.h"

FeaturePaired::FeaturePaired(Feature *a1, Feature *a2) :
		Feature(-1, a1->getName() + " x " + a2->getName(), false), f1(a1), f2(
				a2) {
}

/**
 * Attention: A shallow copy
 */
FeaturePaired::FeaturePaired(const FeaturePaired &o) :
		Feature(o), f1(o.f1), f2(o.f2) {
}

FeaturePaired::~FeaturePaired() {

}

const int FeaturePaired::getFunIndex() const {
	std::cerr << "FeaturePaired error: " << std::endl;
	std::cerr << "A paired feature does not have funIndex.";
	std::cerr << std::endl;
	throw std::exception();
}

int FeaturePaired::getNumOfComp() {
	return 2;
}

int FeaturePaired::getCompOneIndex() {
	return f1->getTableIndex();
}

int FeaturePaired::getCompTwoIndex() {
	return f2->getTableIndex();
}

void FeaturePaired::setIsSelected() {
	isSelected = true;
	f1->setIsNeeded();
	f2->setIsNeeded();
}

void FeaturePaired::setIsNeeded() {
	std::cerr << "FeaturePaired error: " << std::endl;
	std::cerr << "setIsNeeded() is not supported.";
	std::cerr << std::endl;
	throw std::exception();
}

void FeaturePaired::setCompOne(Feature *f) {
	f1 = f;
}

void FeaturePaired::setCompTwo(Feature *f) {
	f2 = f;
}
