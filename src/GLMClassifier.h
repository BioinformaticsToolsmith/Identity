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
 * GLMClassifier.h
 *
 *  Created on: Mar 31, 2020
 *      Author: Dr. Hani Z. Girgis
 */

#ifndef GLMCLASSIFIER_H_
#define GLMCLASSIFIER_H_

#include <iostream>
#include <vector>

#include "ITransformer.h"
#include "Normalizer.h"
#include "SimConverter.h"
#include "FeatureExpander.h"
#include "GLM.h"
#include "BestFirst.h"
#include "Evaluator.h"
#include "StatisticInfo.h"
#include "Feature.h"
#include "GLMPredictor.h"

class GLMClassifier: public ITransformer {
private:
	/**
	 * These are training or validation metrics
	 * They are due to the validation step if it is done
	 */
	double acc = 0.0;
	double sens = 0.0;
	double spec = 0.0;

protected:
	const Matrix *fModelTable;
	const Matrix *lModelTable;

	Matrix *fTrainTable;
	Matrix *fValidateTable;
	Matrix *lTrainTable;
	Matrix *lValidateTable;
	double threshold;
	int threadNum;
	int minFeat;

	std::vector<ITransformer*> *pipe;
	// A list of initial features
	std::vector<Feature*> *fList;
	// A list of optimized features
	std::vector<Feature*> f5;
	// List containing indexes of needed single features
	std::vector<int> indexList;

	// These methods should be overridden in the regressor
	virtual pair<Matrix, std::vector<Feature*> > selectFeatures(Matrix&,
			std::vector<Feature*>&);
	virtual pair<Matrix, GLM*> trainGLM(Matrix&);

	virtual void prepareData();

	virtual void train();
	virtual void validate();

	void clean(std::vector<Feature*>&);

public:
	GLMClassifier(const Matrix*, const Matrix*, double, int, int);
	virtual ~GLMClassifier();
	virtual void start();

	// Methods from ITransformer
	virtual Matrix transform(const Matrix&);
	virtual std::vector<Feature*> getFeatureList();

	virtual void evaluate(Matrix &o, Matrix &p);

	double getAcc() const;
	double getSens() const;
	double getSpec() const;
};

#endif /* GLMCLASSIFIER_H_ */
