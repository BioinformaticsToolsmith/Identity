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
 * ReaderAlignerCoordinator.h
 *
 *  Created on: Nov 15, 2019
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef READERALIGNERCOORDINATOR_H_
#define READERALIGNERCOORDINATOR_H_

#include <vector>
#include <future>

#include "FastaReader.h"
#include "Aligner.h"
#include "Parameters.h"
#include "GLMClassifier.h"
#include "GLMRegressor.h"
#include "SynDataGenerator.h"
#include "AlignerParallel.h"

using namespace std;

class ReaderAlignerCoordinator {
private:
	int workerNum;
	int blockSize;
	char mode;
	double threshold = 0.0;
	bool canRelax;

	void alignFileVsFile1(string, string, string, string, bool);
	void alignFileVsFile2(string, string, string, string, bool);

	template<class V>
	void helper(string, string, string, string, bool, DataGenerator*);

public:
	static const char C = 'c';
	static const char R = 'r';

	ReaderAlignerCoordinator(int, int, char, double, bool);
	virtual ~ReaderAlignerCoordinator();

	void alignAllVsAll(string, string, string);
	void alignQueryVsAll(string, string, string, string);
};

#endif /* READERALIGNERCOORDINATOR_H_ */
