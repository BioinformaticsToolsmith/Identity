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
 * Trainer.cpp
 *
 *  Created on: Feb 25, 2020
 *      Author: Hani Z. Girgis, PhD
 */
#include "SynDataGenerator.h"

SynDataGenerator::SynDataGenerator(string fileName, double t, int threadNumIn,
		std::vector<int> funIndexList) :
		DataGenerator(fileName) {
	threshold = t;
	this->funIndexList = funIndexList;
	threadNum = threadNumIn;
	fillCompositionList();
	generateData();
}

SynDataGenerator::SynDataGenerator(string dbName, string qryName, double t,
		int threadNumIn, std::vector<int> funIndexList) :
		DataGenerator(dbName, qryName, t) {

	threshold = t;
	this->funIndexList = funIndexList;
	threadNum = threadNumIn;
	fillCompositionList();
	generateData();
}

SynDataGenerator::~SynDataGenerator() {
	if (fTable != nullptr) {
		delete fTable;
	}

	if (lTable != nullptr) {
		delete lTable;
	}

	delete[] compositionList;
}

/**
 * Delete data
 * Must be called when all training is finished
 */
void SynDataGenerator::clearData() {
	delete fTable;
	delete lTable;
	fTable = nullptr;
	lTable = nullptr;
}

void SynDataGenerator::fillCompositionList() {
	compositionList = new double[4];
	for (int i = 0; i < 4; i++) {
		compositionList[i] = 0.0;
	}

	int seqNum = block->size();

	// Calculate nucleotide composition vector
	for (int h = 0; h < seqNum; h++) {
		string seq = *(block->at(h).second);
		int len = seq.size();
		for (int i = 0; i < len; i++) {
			switch (seq[i]) {
			case 'A':
				compositionList[0]++;
				break;
			case 'C':
				compositionList[1]++;
				break;
			case 'G':
				compositionList[2]++;
				break;
			case 'T':
				compositionList[3]++;
				break;
			}
		}
	}
	double total = 0.0;
	for (int i = 0; i < 4; i++) {
		total += compositionList[i];
	}
	for (int i = 0; i < 4; i++) {
		compositionList[i] /= total;
	}
}

void SynDataGenerator::generateData() {
	// Determine histogram data type
	if (maxLength <= std::numeric_limits<int8_t>::max()) {
		std::cout << "A histogram entry is 8 bits." << std::endl;
		generateDataHelper<int8_t>();
	} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
		std::cout << "A histogram entry is 16 bits." << std::endl;
		generateDataHelper<int16_t>();
	} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
		std::cout << "A histogram entry is 32 bits." << std::endl;
		generateDataHelper<int32_t>();
	} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
		std::cout << "A histogram entry is 64 bits." << std::endl;
		generateDataHelper<int64_t>();
	} else {
		std::cout << "Trainer warning: Overflow is possible however unlikely.";
		std::cout << std::endl;
		std::cout << "A histogram entry is 64 bits." << std::endl;
		generateDataHelper<int64_t>();
	}
}

/**
 * Build data
 */
template<class V>
void SynDataGenerator::generateDataHelper() {
	std::cout << "Generating data." << std::endl;

	bool canGenerateNegatives = !Util::isEqual(threshold, 0.0);

	// Calculate mutation rates
	double m = Parameters::getMinId();

	// A list containing positive mutation rates
	vector<double> pstvRateList;
	for (double i = 0.99; i >= threshold; i -= 0.01) {
		pstvRateList.push_back(1.0 - i);
	}
	pstvRateList.shrink_to_fit();
	int pstvRateSize = pstvRateList.size();

	// A list containing negative mutation rates
	vector<double> ngtvRateList;
	int ngtvRateSize;
	if (canGenerateNegatives) {
		for (double j = threshold - 0.01; j >= m; j -= 0.01) {
			ngtvRateList.push_back(1.0 - j);
		}
		ngtvRateList.shrink_to_fit();
		ngtvRateSize = ngtvRateList.size();
	}

	int actual = block->size();
	int desired = Parameters::getBlockSize();
	int copyNum = Parameters::getMutPerTemp() / 2.0;
	if (actual < desired) {
		copyNum = copyNum * desired / (double) actual;
	}

	// Generate mutated sequences from each sequence
	KmerHistogram<uint64_t, V> kTable(k);
	KmerHistogram<uint64_t, uint64_t> monoTable(1);
	uint8_t keyList[histogramSize * k];
	kTable.getKeysDigitFormat(keyList);
	const int statNum =
			funIndexList.empty() ?
					StatisticInfo::getInstance()->getStatNum() :
					funIndexList.size();

	// Initialize Matrixes
	int rowNum = 2.0 * copyNum * actual;
	fTable = new Matrix(rowNum, statNum);
	lTable = new Matrix(rowNum, 1);

	if (!canGenerateNegatives) {
		copyNum *= 2;
	}

	int minBlockSize = Parameters::getMutMinBlockSize();
	int maxBlockSize = Parameters::getMutMaxBlockSize();

	bool isSingle = Parameters::getMutSingle();
	bool isBlock = Parameters::getMutBlock();
	bool isTranslocation = Parameters::getMutTranslocation();
	bool isInversion = Parameters::getMutInverstion();

#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < actual; i++) {
		// Set up a mutator
		Mutator mutator(block->at(i).second, maxBlockSize, i, minBlockSize);
		if (isSingle) {
			mutator.enableSinglePoint();
		}
		if (isBlock) {
			mutator.enableBlock();
		}
		if (isTranslocation) {
			mutator.enableTranslocation();
		}
		if (isInversion) {
			mutator.enableInverstion();
		}

		V *h1 = kTable.build(block->at(i).second);
		uint64_t *mono1 = monoTable.build(block->at(i).second);

		// Iterate over different mutation rates
		// Balance around threshold
		// Generate positive examples (identity score above the threshold)
		for (int j = 0; j < copyNum; j++) {
			auto pPstv = mutator.mutateSequence(
					pstvRateList[(i * copyNum + j) % pstvRateSize]);

			V *h2 = kTable.build(pPstv.first);
			uint64_t *mono2 = monoTable.build(pPstv.first);

			Statistician<V> s(histogramSize, k, h1, h2, mono1, mono2,
					compositionList, keyList);
			vector<double> statList;
			statList.reserve(statNum);
			if (funIndexList.empty()) {
				s.calculateAll(statList);
			} else {
				s.calculate(funIndexList, statList);
			}
			int r = canGenerateNegatives ?
					2 * i * copyNum + j : i * copyNum + j;

			fTable->setRow(r, statList);
			lTable->at(r, 0) = pPstv.second;

			delete pPstv.first;
			delete[] h2;
			delete[] mono2;
		}
		// Generate negative examples (identity score below the threshold)
		if (canGenerateNegatives) {
			for (int j = 0; j < copyNum; j++) {
				auto pNgtv = mutator.mutateSequence(
						ngtvRateList[(i * copyNum + j) % ngtvRateSize]);
				V *h2 = kTable.build(pNgtv.first);
				uint64_t *mono2 = monoTable.build(pNgtv.first);

				Statistician<V> s(histogramSize, k, h1, h2, mono1, mono2,
						compositionList, keyList);
				vector<double> statList;
				statList.reserve(statNum);
				if (funIndexList.empty()) {
					s.calculateAll(statList);
				} else {
					s.calculate(funIndexList, statList);
				}
				int r = 2 * i * copyNum + j + copyNum;
				fTable->setRow(r, statList);
				lTable->at(r, 0) = pNgtv.second;

				delete pNgtv.first;
				delete[] h2;
				delete[] mono2;
			}
		}
		delete[] h1;
		delete[] mono1;
	}
}
