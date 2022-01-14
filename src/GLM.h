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
 * glm.h
 *
 * Created on: May 29, 2017
 * Author: Robert Geraghty, The Bioinformatics Toolsmith Laboratory
 * Author: Hani Z. Girgis, PhD,  The Bioinformatics Toolsmith Laboratory
 */

#ifndef SRC_MATRIX_GLM_H_
#define SRC_MATRIX_GLM_H_

#include "Matrix.h"
#include "ITransformer.h"

using namespace std;

typedef double (*output)(double);

class GLM: public ITransformer {
private:
	Matrix weights;
	output o;

public:
	GLM(const Matrix &features, const Matrix &labels, output fun);
	GLM(GLM &other);
	GLM& operator=(const GLM &other);

	virtual Matrix transform(const Matrix &features);

	Matrix getWeights() const;

	virtual std::vector<Feature*> getFeatureList();

	static GLM classifierFactory(const Matrix &f, const Matrix &l);
	static GLM regressorFactory(const Matrix &f, const Matrix &l);

	static GLM* classifierFactoryHeap(const Matrix &f, const Matrix &l);
	static GLM* regressorFactoryHeap(const Matrix &f, const Matrix &l);
};

#endif /* SRC_MATRIX_GLM_H_ */
