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
#include "KmerHistogram.h"
#include "LockFreeQueue.h"
#include "IdentityCalculator.h"

template<class V>
class Aligner {
private:
	IdentityCalculator<V> &identity;
	LockFreeQueue<pair<Block*, bool>, 1000> buffer; // It was 500
	bool canStop = false;
	Block *blockA;
	std::string dlm;
	// If enabled aligner does not align two sequences if they
	// can achieve the minimum identity score.
	bool canReportAll;
	double threshold;
	// Used for relaxing the final filter as threshold - error
	double error = 0.0;
	// Collect all results here
	stringstream *ssPtr;
	// If true write out the content of ssPtr
	bool canWrite = false;

	int k;


//	template<class V>
//	void processBlockHelper();
//	template<class V>
//	void alignSeqVsBlock(std::string *info1, std::string *seq1, Block *blockB,
//			KmerHistogram<uint64_t, V> &kTable,
//			KmerHistogram<uint64_t, uint64_t> &monoTable, int init);
public:
	Aligner(IdentityCalculator<V>&, Block*, string, bool, double, bool);
	virtual ~Aligner();
	void enqueueBlock(pair<Block*, bool>);
	pair<bool, stringstream*> start();
	void stop();
	int getQueueSize();
	void processBlock();
	pair<bool, stringstream*> getResults();
};

#include "Aligner.cpp"

#endif /* ALIGNER_H_ */
