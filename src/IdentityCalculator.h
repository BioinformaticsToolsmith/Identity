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

#ifndef SRC_IDENTITYCALCULATOR_H_
#define SRC_IDENTITYCALCULATOR_H_

#include <vector>
#include <tuple>
#include <algorithm>

#include "DataGenerator.h"
#include "SynDataGenerator.h"
#include "ITransformer.h"
#include "GLMRegressor.h"
#include "GLMPredictor.h"
#include "Parameters.h"
#include "Feature.h"
#include "Matrix.h"
#include "KmerHistogram.h"
#include "Serializer.h"
#include "Util.h"

template<class V>
class IdentityCalculator {
private:
	double threshold;
//	double relaxedThreshold;
	double absError;

//	DataGenerator *g;

	GLMPredictor p;
	int featNum;
	int singleFeatNum;
	std::vector<int> funIndexList;
	int *funIndexArray;

	double *compositionList;
	int k;

	uint8_t *keyList;

	KmerHistogram<uint64_t, V> *kTable;
	KmerHistogram<uint64_t, uint64_t> *monoTable;

	bool canSkip;
	bool canRelax;

protected:
	int kHistSize;
	int monoHistSize;

public:
	IdentityCalculator(DataGenerator*, int, double, bool, bool,
			std::string modelFile = "");
	IdentityCalculator(Serializer&, double, bool, bool);

	virtual ~IdentityCalculator();
	double getError() const;
	int getK() const;

	void freeBlock(std::tuple<V**, uint64_t**, std::string**, int*>, int, int);

	/**
	 * Returns true if it is impossible to obtain the desired identity
	 * score according to the lengths of the two sequences
	 */
	inline bool isImpossible(double l1, double l2, double t) {
//		if (l1 < l2) {
//			std::swap<double>(l1, l2);
//		}
//		bool res = false;
//		if ((l2 / l1) < t) {
//			res = true;
//		}
//		return res;

		double ratio = l1 < l2 ? l1 / l2 : l2 / l1;
		return ratio < t;
	}

	inline double calcRatio(double l1, double l2) {
		return l1 < l2 ? l1 / l2 : l2 / l1;
	}

	/**
	 * One vs. one
	 */
	inline virtual double score(V *kHist1, V *kHist2, uint64_t *monoHist1,
			uint64_t *monoHist2, double ratio, int l1, int l2) {
		// Calculate statistics
		Statistician<V> s(kHistSize, k, kHist1, kHist2, monoHist1, monoHist2,
				compositionList, keyList);
		double res;
		if (canSkip && s.identityMinimum(l1,l2) < threshold) {
			//cout << "Skipping according to filter." << endl;
			res = 0.0;
		} else {
			double data[featNum];
			s.calculate(funIndexArray, singleFeatNum, data);
			// Calculate identity score
			res = p.calculateIdentity(data);

			// In case of error, correct it
			// An identity score cannot be greater than the length ratio
			if (res > ratio) {
				res = ratio;
			}

			// Trim score
			if ((canSkip && res < threshold) || res < 0.0) {
				res = 0.0;
			}
		}
		return res;
	}

	/**
	 * Score one versus many
	 */
	double* score(V *kHist1, V **kHist2List, uint64_t *monoHist1,
			uint64_t **monoHist2List, int listSize, int threadNum, int len1,
			int *len2List);

	/**
	 * Score all versus all
	 * Too-short or too-long pairs may be skipped without evaluation
	 * The threshold may be relaxed according to the regression error
	 */
	Matrix score(V **kHistList, uint64_t **monoHistList, int listSize,
			int threadNum, int *lenList);

// Same as the above method but the matrix is placed on the heap
//	Matrix* scoreHeap(V **kHistList, uint64_t **monoHistList, int listSize,
//			int threadNum, int *lenList);

	std::tuple<V**, uint64_t**, std::string**, int*> unpackBlock(Block *b,
			int threadNum);
	int getKHistSize() const;
	int getMonoHistSize() const;

	bool getCanSkip() const;
	void setCanSkip(bool);
};

#include "IdentityCalculator.cpp"

#endif /* SRC_IDENTITYCALCULATOR_H_ */
