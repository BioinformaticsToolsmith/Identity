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
 * Aligner.cpp
 *
 *  Created on: Nov 15, 2019
 *      Author: Hani Zakaria Girgis, PhD
 *  An instance of this class calculates the All-vs-All on a block of sequences.
 *  It is thread safe designed for parallel execution.
 *
 */

#include "Aligner.h"

/**
 * Block a will NOT be deleted here because it is
 * being processed by other threads as well.
 */
Aligner::Aligner(DataGenerator *d, ITransformer *t, Block *a, string dlmIn,
		bool filter, double cutoff, double e, uint8_t *keyListIn) {
	blockA = a;
	dlm = dlmIn;
	id = t;
	threshold = cutoff;
	error = e;
	keyList = keyListIn;

	histogramSize = d->getHistogramSize();
	k = d->getK();
	maxLength = d->getMaxLength();
	compositionList = d->getCompositionList();

	int alphaSize = Parameters::getAlphabetSize();

	auto featList = t->getFeatureList();
	featNum = featList.size() - 1; // The bias has not been removed yet.
	for (auto f : featList) {
		if (f->getNumOfComp() == 0 && f->getName().compare("constant") != 0) {
			funIndexList.push_back(f->getFunIndex());
		}
	}
	singleFeatNum = funIndexList.size();
	funIndexArray = funIndexList.data();
	// This predictor removes the bias from the feature list
	predictor = GLMPredictor(featList, false);
	// featList is not needed beyond this point
	for (auto f : featList) {
		delete f;
	}
	featList.clear();

	canReportAll = filter;

	ssPtr = new stringstream();
}

Aligner::~Aligner() {
	delete ssPtr;

	if (buffer.size() > 0) {
		std::cerr << "Aligner error: Queue must be empty. " << std::endl;
		std::cerr << "Queue size is: " << buffer.size() << std::endl;
	}
}

pair<bool, stringstream*> Aligner::start() {
	// Keep processing blocks as they are enqueued.
	while (true) {
		if (buffer.size() > 0) {
			processBlock();
		} else if (canStop && buffer.size() == 0) {
			break;
		}
	}

	return std::make_pair(canWrite, ssPtr);
}

/**
 * Thread safe
 * Note: This block and its contents will be deleted
 * after processing it.
 */
void Aligner::enqueueBlock(pair<Block*, bool> p) {
	buffer.push(p);
}

void Aligner::processBlock() {
	// Determine histogram data type
	if (maxLength <= std::numeric_limits<int8_t>::max()) {
		processBlockHelper<int8_t>();
	} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
		processBlockHelper<int16_t>();
	} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
		processBlockHelper<int32_t>();
	} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
		processBlockHelper<int64_t>();
	} else {
		std::cout << "Aligner warning: Overflow is possible however unlikely.";
		std::cout << std::endl;
		std::cout << "A histogram entry is 64 bits." << std::endl;
		processBlockHelper<int64_t>();
	}
}

/**
 * Thread safe
 */
template<class V>
void Aligner::processBlockHelper() {
	// Get the front of the the queue, but do not pop it yet.
	auto blockB = buffer.front().first;
	int sizeA = blockA->size();
	int sizeB = blockB->size();

	static KmerHistogram<uint64_t, V> kTable(k);
	static KmerHistogram<uint64_t, uint64_t> monoTable(1);

	for (int j = 0; j < sizeA; j++) {
		int init = 0;
		// If the two blocks have the same contents.
		if (buffer.front().second) {
			init = j + 1;
		}
		auto p1 = blockA->at(j);
		alignSeqVsBlock(p1.first, p1.second, blockB, kTable, monoTable, init);
	}

	// Pop the block and free its memory
	FastaReader::deleteBlock(blockB);
	buffer.pop();
}

template<class V>
void Aligner::alignSeqVsBlock(std::string *info1, std::string *seq1,
		Block *blockB, KmerHistogram<uint64_t, V> &kTable,
		KmerHistogram<uint64_t, uint64_t> &monoTable, int init) {

	V *h1 = kTable.build(seq1);
	uint64_t *mono1 = monoTable.build(seq1);

	int l1 = seq1->size();
	int sizeB = blockB->size();
	double data[featNum];

	for (int hani = init; hani < sizeB; hani++) {
		auto p2 = blockB->at(hani);
		string *seq2 = p2.second;
		int l2 = seq2->size();

		if (!canReportAll
				&& (std::min(l1, l2) / (double) std::max(l1, l2) < threshold)) {
			continue;
		}

		V *h2 = kTable.build(seq2);
		uint64_t *mono2 = monoTable.build(seq2);

		Statistician<V> s(histogramSize, k, h1, h2, mono1, mono2,
				compositionList, keyList);

		s.calculate(funIndexArray, singleFeatNum, data);

		double res = predictor.calculateIdentity(data);

		if (canReportAll || res >= threshold - error) {
			canWrite = true;

			if (res > 1.0) {
				res = 1.0;
			} else if (res < 0.0) {
				res = 0.0;
			}

			(*ssPtr) << *info1 << dlm << *p2.first << dlm
					<< std::setprecision(4) << res << std::endl;
		}

		delete[] h2;
		delete[] mono2;
	}

	delete[] h1;
	delete[] mono1;
}

int Aligner::getQueueSize() {
	return buffer.size();
}

/**
 * Thread safe
 */
void Aligner::stop() {
	canStop = true;
}
