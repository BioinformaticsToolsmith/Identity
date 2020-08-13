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
 * ITransformer.cpp
 *
 *  Created on: May 20, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#include "ITransformer.h"

/**
 * Make a deep copy of each feature in l
 * It sets component pointers according to the their locations
 * in the new vector
 */
std::vector<Feature*> ITransformer::copy(const std::vector<Feature*> &l) {
	int s = l.size();
	// Make sure that l is sorted as singles, squared, paired
	for (int i = 0; i < s - 1; i++) {
		if (l[i]->getNumOfComp() > l[i + 1]->getNumOfComp()) {
			std::cerr << "ITransformer error: " << std::endl;
			std::cerr
					<< "Feature list must be sorted (singles, squared, paired).";
			std::cerr << std::endl;
			throw std::exception();
		}
	}

	// Make sure every feature has an index in the table
	for (int i = 0; i < s; i++) {
		if (l[i]->getTableIndex() == -1) {
			std::cerr << "ITransformer error: " << std::endl;
			std::cerr << l[i]->getName() << " has uninitialized table index.";
			std::cerr << std::endl;
			throw std::exception();
		}
	}

	// Perform deep copy
	std::vector<Feature*> r(s, nullptr);
	for (int i = 0; i < s; i++) {
		auto f = l[i];
		int c = f->getNumOfComp();
		if (c == 2) {
			auto t = new FeaturePaired(*static_cast<FeaturePaired*>(f));
			t->setCompOne(r.at(f->getCompOneIndex()));
			t->setCompTwo(r.at(f->getCompTwoIndex()));
			r[i] = t;
		} else if (c == 1) {
			auto t = new FeatureSquared(*static_cast<FeatureSquared*>(f));
			t->setCompOne(r.at(f->getCompOneIndex()));
			r[i] = t;
		} else if (c == 0) {
			r[i] = new Feature(*l[i]);
		} else {
			std::cerr << "ITransformer error: " << std::endl;
			std::cerr << "Unexpected component number: " << c;
			std::cerr << std::endl;
			throw std::exception();
		}
	}
	return r;
}
