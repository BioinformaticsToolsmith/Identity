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
 * Aligner.h
 *
 *  Created on: Nov 15, 2019
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <iostream>
#include <sstream>
#include <vector>
#include "Util.h"
#include "FastaReader.h"
#include "ITransformer.h"
#include "DataGenerator.h"
#include "KmerHistogram.h"
#include "Statistician.h"
#include "Feature.h"
#include "LockFreeQueue.h"
#include "GLMPredictor.h"

class Aligner {
private:
	LockFreeQueue<pair<Block*, bool>, 500> buffer;
	bool canStop = false;

	Block *blockA;
	std::string dlm;
	ITransformer *id;

	// Needed for histograms and statisticians
	// Get them from data generator
	int histogramSize;
	int k;
	int64_t maxLength;
	double *compositionList;

	std::vector<int> funIndexList;
	int *funIndexArray;
	int singleFeatNum;
	int featNum;

	// If enabled aligner does not align two sequences if they
	// can achieve the minimum identity score.
	bool isLengthFilter;

	double threshold;

	// Used for relaxing the final filter as threshold - error
	double error = 0.0;

	// Collect all results here
	stringstream *ssPtr;
	// If true write out the content of ssPtr
	bool canWrite = false;

	GLMPredictor predictor;

	uint8_t *keyList;

	void processBlock();
	template<class V>
	void processBlockHelper();
	template<class V>
	void alignSeqVsBlock(std::string *info1, std::string *seq1, Block *blockB,
			KmerHistogram<uint64_t, V> &kTable,
			KmerHistogram<uint64_t, uint64_t> &monoTable, /*uint8_t *keyList,*/
			int init);

public:
	Aligner(DataGenerator*, ITransformer*, Block*, string, bool, double,
			double e, uint8_t*);
	virtual ~Aligner();
	void enqueueBlock(pair<Block*, bool>);
	pair<bool, stringstream*> start();
	void stop();
	int getQueueSize();
};

#endif /* ALIGNER_H_ */
