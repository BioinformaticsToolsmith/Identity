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
 * GLMClassifier.h
 *
 *  Created on: Mar 31, 2020
 *      Author: Dr. Hani Z. Girgis
 */

#ifndef GLMREGRESSOR_H_
#define GLMREGRESSOR_H_

#include <iostream>
#include <vector>
#include <limits>

#include "GLM.h"
#include "BestFirst.h"
#include "Evaluator.h"
#include "GLMClassifier.h"
#include "GLMPredictor.h"

class GLMRegressor: public GLMClassifier {
private:
	/**
	 * These are training or validation metrics
	 * They are due to the validation step if it is done
	 */
	double absError = 0.0;
	// double rlxError = 0.0;
	double sqrError = 0.0;

protected:
	virtual void prepareData();

	virtual pair<Matrix, std::vector<Feature*> > selectFeatures(Matrix&,
			std::vector<Feature*>&);

	virtual pair<Matrix, GLM*> trainGLM(Matrix&);

public:
	GLMRegressor(const Matrix*, const Matrix*, double, int, int);
	virtual ~GLMRegressor();

	virtual void evaluate(Matrix &o, Matrix &p);
	double getAbsError() const;
	double getSqrError() const;
};

#endif /* GLMCLASSIFIER_H_ */
