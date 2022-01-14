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
 * Util.h
 *
 *  Created on: Apr 18, 2020
 *      Author: Dr. Hani Z. Girgis
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <limits>
#include <math.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Parameters.h"

class Util {
public:
	static inline bool isEqual(double d1, double d2) {
		bool r = false;
		if (fabs(d1 - d2) < std::numeric_limits<double>::epsilon()) {
			r = true;
		}
		return r;
	}

	static inline uint8_t* makeKeyList(int histSize, int k) {
		uint8_t *keyList = new uint8_t[histSize * k];
		int alphaSize = Parameters::getAlphabetSize();
		for (int c = k - 1; c >= 0; c--) {
			for (int r = 0; r < histSize; r++) {
				keyList[(r * k) + k - 1 - c] = (r
						/ ((uint64_t) pow(alphaSize, c))) % alphaSize;
			}
		}
		return keyList;
	}

	static inline bool doesFileExist(std::string fileName) {
		std::ifstream in(fileName);
		bool r = false;
		if (in.good()) {
			r = true;
		}
		return r;
	}

	static double calculateMedian(std::vector<double> v) {
		std::sort(v.begin(), v.end(), [](double i, double j) -> bool {
			return i > j;
		});
		return v.at(v.size() / 2);
	}

	static double calculateMean(std::vector<double> &v) {
		double mean = 0.0;
		for (double d : v) {
			mean += d;
		}
		return mean / v.size();
	}

	static double calculateSTDev(std::vector<double> &v, double mean) {
		double sum = 0.0;
		for (double d : v) {
			double diff = d - mean;
			sum += diff * diff;
		}
		return sqrt(sum / v.size());
	}

	template<class V>
	static bool isAllZeros(V *hist, int histSize) {
		bool r = true;
		for (int i = 0; i < histSize; i++) {
			if (hist[i] > 0.0) {
				r = false;
				break;
			}
		}
		return r;
	}
};

#endif /* UTIL_H_ */
