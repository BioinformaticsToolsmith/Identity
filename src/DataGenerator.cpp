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
 * DataGenerator.cpp
 *
 *  Created on: Apr 23, 2020
 *      Author: Dr. Hani Zakaria Girgis
 */

#include "DataGenerator.h"

/**
 * This constructor should be used for all vs. all
 */
DataGenerator::DataGenerator(std::string fileName, int blockSize) {
	FastaReader reader(fileName, blockSize);
	block = reader.read();

	calculateK();
	calculateHistSize();

	// Find maximum length
	int seqNum = block->size();
	for (int h = 0; h < seqNum; h++) {
		uint64_t len = block->at(h).second->size();
		if (len > maxLength) {
			maxLength = len;
		}
	}
	// The length was estimated based on a subset not the entire set
	maxLength = 2 * maxLength;
}

/**
 * This constructor should be used for searching for query sequences
 * in a database. It is purpose is to select training sample similar
 * in length to the query sequences.
 */
DataGenerator::DataGenerator(std::string dbName, std::string qryFile,
		double threshold) {

	FastaReader qryReader(qryFile, Parameters::getBlockSize());
	Block *qryBlock = qryReader.read();

	FastaReader dbReader(dbName, Parameters::getBlockSize());
	block = dbReader.read();

	// Move query sequence(s) to the database block
	for (auto p : *qryBlock) {
		block->push_back(p);
	}

	if (threshold > 0.0) {

		int min = std::numeric_limits<int>::max();
		int max = -1;

		for (auto p : *qryBlock) {
			int s = p.second->size();

			if (s < min) {
				min = s;
			}

			if (s > max) {
				max = s;
			}
		}
		min = min * threshold;
		max = max / threshold;

		// Filter the database block according to length
		int seqNum = block->size();
		for (int e = seqNum - 1; e >= 0; e--) {
			auto p = block->at(e);
			int s = p.second->size();
			if (s < min || s > max) {
				delete p.first;
				delete p.second;
				block->at(e) = block->back();
				block->pop_back();
			}
		}

		if (block->empty()) {
			std::cerr
					<< "Data generator error: No sequences left after filtering."
					<< std::endl;
			throw std::exception();
		}

		maxLength = max;
	} else {
		// Find maximum length
		int seqNum = block->size();
		for (int h = 0; h < seqNum; h++) {
			uint64_t len = block->at(h).second->size();
			if (len > maxLength) {
				maxLength = len;
			}
		}
		// The length was estimated based on a subset not the entire set
		maxLength = 2 * maxLength;
	}

	calculateK();
	calculateHistSize();

	// Free up memory
	delete qryBlock;
}

void DataGenerator::calculateK() {
	// Calculate average length
	double sum = 0.0;
	for (auto p : *block) {
		sum += p.second->size();
	}

	average = round(sum / block->size());
	std::cout << "Average: " << average << std::endl;
	// Calculate k
	k = ceil(log(average) / log(Parameters::getAlphabetSize()))
			- Parameters::getKRelax();
	if (k < 2) {
		std::cerr << "DataGenerator warning: K is too small. ";
		std::cerr << "K is set to 2." << std::endl;
		k = 2;
	}
	std::cout << "K: " << k << std::endl;
}

void DataGenerator::calculateHistSize() {
	// Calculate histogram size
	histogramSize = pow(4, k);
	std::cout << "Histogram size: " << histogramSize << std::endl;
}

DataGenerator::~DataGenerator() {
	FastaReader::deleteBlock(block);
}

const Matrix* DataGenerator::getFeatures() const {
	return fTable;
}

double* DataGenerator::getCompositionList() const {
	return compositionList;
}

int DataGenerator::getHistogramSize() const {
	return histogramSize;
}

int DataGenerator::getK() const {
	return k;
}

int64_t DataGenerator::getMaxLength() const {
	return maxLength;
}

const Matrix* DataGenerator::getLabels() const {
	return lTable;
}
