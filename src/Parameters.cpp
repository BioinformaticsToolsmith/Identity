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
 * Parameters.cpp
 *
 *  Created on: Nov 30, 2019
 *      Author: Hani Z. Girgis, PhD
 */

#include "Parameters.h"

/**
 * Default parameters
 * These are meant for the internal use by the program not the user.
 */
int Parameters::MODE = DNA;
double Parameters::MIN_ID = 0.0;
int Parameters::MUT_PER_TEMP = 10;
bool Parameters::IS_SINGLE = true;
bool Parameters::IS_BLOCK = true;
bool Parameters::IS_TRANSLOCATION = false;
bool Parameters::IS_INVERSTION = false;
int Parameters::MIN_BLOCK_SIZE = 2;
int Parameters::MAX_BLOCK_SIZE = 5;
int Parameters::BLOCK_SIZE = 1000;
int Parameters::K_RELAX = 1; //2;

// These are not used
int Parameters::MIN_FEAT_NUM = 3;
int Parameters::MAX_FEAT_NUM = 5;
double Parameters::DELTA_R = 0.000025; // Unused
double Parameters::DELTA_C = 0.001; // Unused

double Parameters::TRANSLOCATION_FACTOR = 1.0;
double Parameters::INVERSION_FACTOR = 1.0;

int Parameters::MS_ITR = 100;
int Parameters::MS_BANDWIDTH_BLOCK = 10000; // 1000
int Parameters::MS_BLOCK = 25000;
int Parameters::MS_V_BLOCK = 100000;
int Parameters::MS_READ_MORE = 100000;
int Parameters::MS_PASS_NUM = 10;
double Parameters::MS_BANDWIDTH_THRESHOLD = 0.7;
int Parameters::MS_MAX_MATRIX_SIZE = 46340;
int Parameters::MS_BANDWIDTH_ITERATIONS = 3;
int Parameters::MS_PRINT_BLOCK = 50000;
double Parameters::MS_SLACK_MAX = 0.02;

Parameters::Parameters() {
	checkMutPerTemp();
	checkMinId();
}

Parameters::~Parameters() {
}

double Parameters::getTranslocationFactor() {
	return TRANSLOCATION_FACTOR;
}

double Parameters::getInversionFactor() {
	return INVERSION_FACTOR;
}

int Parameters::getMode() {
	return MODE;
}

bool Parameters::isDNA() {
	bool r = false;
	if (MODE == DNA) {
		r = true;
	}
	return r;
}

void Parameters::setMode(int newMode) {
	MODE = newMode;
}

char Parameters::getUnknown() {
	char unknown;
	if (MODE == DNA) {
		unknown = 'N';
	} else if (MODE == PROTEIN) {
		unknown = 'X';
	} else {
		std::cerr << "Unknown mode: " << MODE << std::endl;
		throw std::exception();
	}
	return unknown;
}

int Parameters::getAlphabetSize() {
	int s;
	if (MODE == DNA) {
		s = 4;
	} else if (MODE == PROTEIN) {
		s = 22;
	} else {
		std::cerr << "Unknown mode: " << MODE << std::endl;
		throw std::exception();
	}
	return s;
}

// Parameter of Trainer
double Parameters::getMinId() {
	return MIN_ID;
}

void Parameters::checkMinId() {
	if (MIN_ID > 1.0) {
		std::cerr << "The minimum identity cannot be greater than 1."
				<< std::endl;
		throw std::exception();
	}

	if (MIN_ID < 0.0) {
		std::cerr << "The minimum identity cannot be less than 0." << std::endl;
		throw std::exception();
	}
}

// Parameter of Trainer
void Parameters::setMinId(double newMin) {
	MIN_ID = newMin;
	checkMinId();
}

// Parameter of Trainer
int Parameters::getMutPerTemp() {
	return MUT_PER_TEMP;
}

void Parameters::checkMutPerTemp() {
	// Number of mutated sequences per template must be even,
	// so we can generate equal number of positive and negative pairs.
	if (MUT_PER_TEMP % 2 != 0) {
		std::cerr << "Number of mutated sequences per template must be even.";
		std::cerr << std::endl;
		throw std::exception();
	}
}

// Parameter of Trainer
void Parameters::setMutPerTemp(int newNum) {
	MUT_PER_TEMP = newNum;
	checkMutPerTemp();
}

// Parameter of Reader
int Parameters::getBlockSize() {
	return BLOCK_SIZE;
}

// Parameter of Reader
void Parameters::setBlockSize(int newSize) {
	BLOCK_SIZE = newSize;
}

// Parameter of Mutator
void Parameters::enableMutSingle() {
	IS_SINGLE = true;
}

// Parameter of Mutator
void Parameters::enableMutBlock() {
	IS_BLOCK = true;
}

// Parameter of Mutator
void Parameters::enableMutTranslocation() {
	IS_TRANSLOCATION = true;
}

// Parameter of Mutator
void Parameters::enableMutInverstion() {
	IS_INVERSTION = true;
}

// Parameter of Mutator
bool Parameters::getMutSingle() {
	return IS_SINGLE;
}

// Parameter of Mutator
bool Parameters::getMutBlock() {
	return IS_BLOCK;
}

// Parameter of Mutator
bool Parameters::getMutTranslocation() {
	return IS_TRANSLOCATION;
}

// Parameter of Mutator
bool Parameters::getMutInverstion() {
	return IS_INVERSTION;
}

// Parameter of Mutator
int Parameters::getMutMinBlockSize() {
	return MIN_BLOCK_SIZE;
}

// Parameter of Mutator
int Parameters::getMutMaxBlockSize() {
	return MAX_BLOCK_SIZE;
}

// Parameter of Mutator
void Parameters::setMutMinBlockSize(int newSize) {
	MIN_BLOCK_SIZE = newSize;
}

// Parameter of Mutator
void Parameters::setMutMaxBlockSize(int newSize) {
	MAX_BLOCK_SIZE = newSize;
}

int Parameters::getMinFeatNum() {
	return MIN_FEAT_NUM;
}

void Parameters::setMinFeatNum(int n) {
	MIN_FEAT_NUM = n;
}

int Parameters::getMaxFeatNum() {
	return MAX_FEAT_NUM;
}

void Parameters::setMaxFeatNum(int n) {
	MAX_FEAT_NUM = n;
}

double Parameters::getDeltaC() {
	return DELTA_C;
}

void Parameters::setDeltaC(double deltaC) {
	if (deltaC < 0.0 || deltaC > 1.0) {
		std::cerr << "A classifier delta must be between 0 & 1. ";
		std::cerr << "Received: " << deltaC << std::endl;
		throw std::exception();
	}
	DELTA_C = deltaC;
}

double Parameters::getDeltaR() {
	return DELTA_R;
}

void Parameters::setDeltaR(double deltaR) {
	if (deltaR < 0.0 || deltaR > 1.0) {
		std::cerr << "A regression delta must be between 0 & 1. ";
		std::cerr << "Received: " << deltaR << std::endl;
		throw std::exception();
	}
	DELTA_R = deltaR;
}

int Parameters::getKRelax() {
	return K_RELAX;
}

void Parameters::setKRelax(int kRelax) {
	K_RELAX = kRelax;
}

int Parameters::getMsItr() {
	return MS_ITR;
}

int Parameters::getMsBandwidthBlock() {
	return MS_BANDWIDTH_BLOCK;
}

double Parameters::getMsBandwidthThreshold() {
	return MS_BANDWIDTH_THRESHOLD;
}

//double Parameters::getMsReliableNumerator() {
//	return MS_RELIABLE_NUMERATOR;
//}

//double Parameters::getMsReliableThreshold() {
//	return MS_RELIABLE_THRESHOLD;
//}

int Parameters::getMsMaxMatrixSize() {
	return MS_MAX_MATRIX_SIZE;
}

int Parameters::getMsBandwidthIterations() {
	return MS_BANDWIDTH_ITERATIONS;
}

int Parameters::getMsBlock() {
	return MS_BLOCK;
}

int Parameters::getMsVBlock() {
	return MS_V_BLOCK;
}

int Parameters::getMsPassNum() {
	return MS_PASS_NUM;
}

double Parameters::getMsSlackMax() {
	return MS_SLACK_MAX;
}

int Parameters::getMsPrintBlock() {
	return MS_PRINT_BLOCK;
}
