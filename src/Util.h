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
 * Util.h
 *
 *  Created on: Apr 18, 2020
 *      Author: Dr. Hani Z. Girgis
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <limits>
#include <math.h>

class Util {
public:
	static inline bool isEqual(double d1, double d2) {
		bool r = false;
		if (fabs(d1 - d2) < std::numeric_limits<double>::epsilon()) {
			r = true;
		}
		return r;
	}
};

#endif /* UTIL_H_ */
