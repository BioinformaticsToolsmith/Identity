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
 * Mutator.h
 *
 *  Created on: Oct 23, 2019
 *      Author: Hani Zakaria Girgia, PhD
 * Notes:
 * I have 8 mutation types
 * 1. Single insertion
 * 2. Single deletion
 * 3. Single change
 * 4. Block insertion
 * 5. Block deletion
 * 6. Duplication
 * 7. Reversion
 * 8. Translocation
 *
 * ToDo: Write the algorithm here!
 *
 */

#ifndef MUTATOR_H_
#define MUTATOR_H_

#include <string>
#include <math.h>
#include <vector>
// #include <stdlib.h> // rand()
#include <random>
#include <algorithm> // count()

#include "Parameters.h"

using namespace std;

class Mutator {
private:
	const string * oSequence;
	vector<int> * mutationList;
	bool ownCompositionList;
	vector<double> * compositionList;
	// A list holding pairs of valid segments, which do not include N or X.
	vector< pair<int,int> > * segmentList;

	// When the block size is small (e.g. < 5 nucleotides), the single point
	// mutation model is very realistic
	int maxBlock = 0;
	int minBlock;
	double aLimit, cLimit, gLimit, tLimit;
	char unknown;
	int effectiveLength = 0;

	default_random_engine * gen;
	uniform_real_distribution<> * zeroOneRand;

	double translocationFactor;
	double inversionFactor;

	char getRandomNucleotide();
	void makeCompositionList();
	void help(int, int , int);

public:
	enum Mutation {INSERTION, DELETION, MISMATCH, B_INSERTION,
		B_DELETION, DUPLICATION, INVERSION, TRANSLOCATION};

	Mutator(const string *, int, int , int minBlockIn = 2);
	Mutator(const string *, int, int , vector<double> *,int minBlockIn = 2);
	virtual ~Mutator();
	void enableSinglePoint();
	void enableBlock();
	void enableInverstion();
	void enableTranslocation();

	pair<string*, double> mutateSequence(double);
};

#endif /* MUTATOR_H_ */
