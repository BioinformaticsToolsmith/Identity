/*
 Identity 2.0 calculates DNA sequence identity scores rapidly without alignment.

 Copyright (C) 2020-2022 Hani Z. Girgis, PhD

 Academic use: Affero General Public License version 1.

 Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 */

/*
 * IdentityCalculator.h
 *
 *  Created on: Dec 24, 2020
 *      Author: Hani Zakaria Girgis, PhD
 *     Purpose: This class is all is needed to calculated identity scores
 */

#ifndef SRC_IDENTITYCALCULATOR1_H_
#define SRC_IDENTITYCALCULATOR1_H_

#include "IdentityCalculator.h"

template<class V>
class IdentityCalculator1: public IdentityCalculator<V> {
public:
	IdentityCalculator1(DataGenerator*, int, double, bool, bool);
	virtual ~IdentityCalculator1();

	/**
	 * One vs. one: Exact and identical match
	 */
	inline virtual double score(V *kHist1, V *kHist2, uint64_t *monoHist1,
			uint64_t *monoHist2) {
		bool isIdentical = true;
		for (int i = 0; i < IdentityCalculator<V>::kHistSize; i++) {
			if (kHist1[i] != kHist2[i]) {
				isIdentical = false;
				break;
			}
		}

		return isIdentical ? 1.0 : 0.0;
	}
};

#include "../src/IdentityCalculator1.cpp"

#endif /* SRC_IDENTITYCALCULATOR1_H_ */
