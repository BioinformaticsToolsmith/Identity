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
 * GLMClassifier.cpp
 *
 *  Created on: Mar 31, 2020
 *      Author: Dr. Hani Z. Girgis
 *
 */
#include "GLMRegressor.h"

GLMRegressor::GLMRegressor(const Matrix *f, const Matrix *l, double t, int c,
		int m) :
		GLMClassifier(f, l, t, c, m) {
}

GLMRegressor::~GLMRegressor() {

}

/**
 * Fill the four matrices here of positive examples
 */
void GLMRegressor::prepareData() {
	std::cout << "Preparing data ..." << std::endl;

	// Determine matrices sizes
	int numPstv = 0;
	int numRow = lModelTable->getNumRow();

	for (int i = 0; i < numRow; i++) {
		if (lModelTable->item(i, 0) >= threshold) {
			numPstv++;
		}
	}

	// Make sure we have positive  examples
	if (numPstv == 0) {
		std::cerr << "GLMRegressor error: No positives." << std::endl;
		throw std::exception();
	}

	// numPstv needs to be even
	if (numPstv % 2 != 0) {
		numPstv--;
	}

	// Build two balanced sets for training and validation
	int p = 0; // Number of collected training and validation positive examples
	int t = 0; // Number of collected training examples
	int v = 0; // Number of collected validation examples
	int i = 0; // Index in model table

	int s = numPstv / 2.0;
	int trainList[s]; // Array holding training indexes
	int validateList[s]; // Array holding validation indexes

	while (p < numPstv && i < numRow) {
		double id = lModelTable->at(i, 0);
		if (id >= threshold) {
			if (p % 2 == 0) {
				trainList[t] = i;
				t++;
			} else {
				validateList[v] = i;
				v++;
			}
			p++;
		}
		i++;
	}

	// Feature and label matrices
	fTrainTable = fModelTable->subMatrix(trainList, s);
	lTrainTable = lModelTable->subMatrix(trainList, s);
	fValidateTable = fModelTable->subMatrix(validateList, s);
	lValidateTable = lModelTable->subMatrix(validateList, s);

	std::cout << "\tPositive examples: " << p << std::endl;
	std::cout << "\tTraining size: " << t << std::endl;
	std::cout << "\tValidation size: " << v << std::endl;
}

pair<Matrix, std::vector<Feature*> > GLMRegressor::selectFeatures(Matrix &t4,
		std::vector<Feature*> &f4) {

	auto isNewBetter = [](double newV, double oldV) {
		return (oldV - newV > 0.000025) ? true : false; // Original
	};

	BestFirst<GLM> selector(t4, *lTrainTable, f4, GLM::regressorFactory,
			Evaluator::mse, isNewBetter, false, threadNum, minFeat,
			std::numeric_limits<double>::infinity());

	return std::make_pair(selector.transform(t4), selector.getFeatureList());
}

pair<Matrix, GLM*> GLMRegressor::trainGLM(Matrix &t5) {
	GLM *glm = GLM::regressorFactoryHeap(t5, *lTrainTable);
	Matrix t6 = glm->transform(t5);
	return std::make_pair(t6, glm);
}

double GLMRegressor::getAbsError() const {
	return absError;
}

double GLMRegressor::getSqrError() const {
	return sqrError;
}

void GLMRegressor::evaluate(Matrix &o, Matrix &p) {
	absError = Evaluator::mae(o, p);
	//rlxError = Evaluator::maeForMeshclust(o, p, threshold);
	sqrError = Evaluator::mse(o, p);

	std::cout << "\tMAE: " << absError << std::endl;
	//std::cout << "\tRLX: " << rlxError << std::endl;
	std::cout << "\tMSE: " << sqrError << std::endl;

}
