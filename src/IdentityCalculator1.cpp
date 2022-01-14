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
template<class V>
IdentityCalculator1<V>::IdentityCalculator1(DataGenerator *generator,
		int threadNum, double t, bool skip, bool relax) :
		IdentityCalculator<V>(generator, threadNum, t, skip, relax) {
}

template<class V>
IdentityCalculator1<V>::~IdentityCalculator1() {

}

