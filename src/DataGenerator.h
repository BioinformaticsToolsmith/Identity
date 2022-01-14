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
 * DataGenerator.h
 *
 *  Created on: Apr 23, 2020
 *      Author: Dr. Hani Zakaria Girgis
 */

#ifndef DATAGENERATOR_H_
#define DATAGENERATOR_H_

#include <string>
#include "FastaReader.h"
#include "KmerHistogram.h"
#include "Matrix.h"
#include "Parameters.h"

class DataGenerator {
private:
	void calculateK();
	void calculateHistSize();

protected:
	Block *block;
	int histogramSize;
	int k;
	const int kRelax = 2;
	int64_t maxLength = 0;
	int64_t average = 0;

	// List of nucleotide composition of generated sequences
	double *compositionList;

	// Model of all data (single statistics only)
	Matrix *fTable = nullptr;
	// Model of associated labels (identity scores)
	Matrix *lTable = nullptr;

public:
	DataGenerator(std::string, int blockSize = Parameters::getBlockSize());
	DataGenerator(std::string, std::string, double);
	virtual ~DataGenerator();

	// Free memory used by the two tables before the object is destroyed
	virtual void clearData() = 0;

	// A Subclass implements this method using different
	// methods for generating sequence pairs with their labels
	virtual void generateData() = 0;

	const Matrix* getFeatures() const;
	const Matrix* getLabels() const;
	double* getCompositionList() const;
	int getHistogramSize() const;
	int getK() const;
	int64_t getMaxLength() const;
};

#endif /* DATAGENERATOR_H_ */
