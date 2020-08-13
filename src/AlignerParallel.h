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
 * AlignerParallel.h
 *
 *  Created on: Jul 20, 2020
 *      Author: Hani Z. Girgis, PhD
 *     Purpose: The single instance of this class aligns one block
 *     versus multiple blocks.
 *     There should be one instance of this class
 */

#ifndef SRC_ALIGNERPARALLEL_H_
#define SRC_ALIGNERPARALLEL_H_

#include <string>
#include <iostream>
#include <vector>
#include <omp.h> // To get thread identifier
#include <future>
#include <atomic>
#include <tuple>

#include "KmerHistogram.h"
#include "Statistician.h"
#include "Parameters.h"
#include "ITransformer.h"
#include "GLMPredictor.h"

typedef std::vector<std::pair<std::string*, std::string*> > Block;
typedef std::vector<std::vector<pair<std::string*, double> >*> Result;

template<class V>
class AlignerParallel {
private:
	// Block A Data
	V **kHistList;
	uint64_t **monoHistList;
	std::string **infoList;
	int *lenList;
	int sizeA = 0;

	// Number of threads to used for processing two blocks
	int threadNum;

	int featNum;
	int singleFeatNum;

	bool isLengthFilter;
	double threshold;
	double relaxThreshold;

	KmerHistogram<uint64_t, V> *kTable;
	KmerHistogram<uint64_t, uint64_t> *monoTable;

	double *compositionList;
	int k;
	int histSize;
	uint8_t *keyList;

	std::vector<int> funIndexList;
	int *funIndexArray;
	GLMPredictor predictor;

	std::string dlm;

	bool isInitialized;

	std::ofstream out;

	void clearAMemory(V**, uint64_t**, std::string**, int*, int);
	void output(std::vector<std::vector<pair<std::string*, double> >*>*,
			std::string*);

	Result* makeEmptyResult(int, int);

public:
	AlignerParallel(int, int, double, double, bool, double*, ITransformer*,
			std::string, int, std::string);
	virtual ~AlignerParallel();
	int getThreadNum() const;
	void setThreadNum(int threadNum);
	std::tuple<V**, uint64_t**, std::string**, int*> unpackBlock(Block*);
	void setBlockA(Block*, bool);
	void processBlockB(Block*);
	bool isDone();
};

#include "AlignerParallel.cpp"

#endif /* SRC_ALIGNERPARALLEL_H_ */
