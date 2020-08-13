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
 * glm.cpp
 *
 * Created on: May 29, 2017
 * Modified on: April 18, 2020
 * Author: Robert Geraghty, The Bioinformatics Toolsmith Laboratory
 * Author: Dr. Hani Z. Girgis, The Bioinformatics Toolsmith Laboratory
 */

#include "GLM.h"
#include "Matrix.h"

#include <math.h>
#include <iostream>

/**
 * f: Feature matrix
 * l: Label matrix
 * fun: output function with default for classification between 0 and 1
 */
GLM::GLM(const Matrix &f, const Matrix &l, output fun) {
	o = fun;
	Matrix t = ~f;
	Matrix s = t * f;
	weights = s.pseudoInverse() * t * l;
}

GLM::GLM(GLM &other) {
	weights = other.weights;
	o = other.o;
}

GLM& GLM::operator=(const GLM &other) {
	weights = other.weights;
	o = other.o;

	return *this;
}

Matrix GLM::transform(const Matrix &features) {
	Matrix labels = features * weights;
	// ToDo: Add pragma
	for (int i = 0; i < labels.getNumRow(); i++) {
		labels.at(i, 0) = o(labels.at(i, 0));
	}
	return labels;
}

Matrix GLM::getWeights() const {
	return weights;
}

/**
 * Makes a GLM instance with a binary output function for classification
 */
GLM GLM::classifierFactory(const Matrix &f, const Matrix &l) {
	return GLM(f, l, [](double x) {
		return x >= 0.5 ? 1.0 : 0.0;
	});
}

/**
 * Makes a GLM instance with a linear output function for regression
 */
GLM GLM::regressorFactory(const Matrix &f, const Matrix &l) {
	return GLM(f, l, [](double x) {
		return x;
	});
}

/**
 * Makes a GLM instance with a binary output function for classification
 * This instance is allocated on the heap. It is the clients responsibility
 * to free its memory.
 */
GLM* GLM::classifierFactoryHeap(const Matrix &f, const Matrix &l) {
	return new GLM(f, l, [](double x) {
		return x >= 0.5 ? 1.0 : 0.0;
	});
}

/**
 * Makes a GLM instance with a linear output function for regression
 * This instance is allocated on the heap. It is the clients responsibility
 * to free its memory.
 */
GLM* GLM::regressorFactoryHeap(const Matrix &f, const Matrix &l) {
	return new GLM(f, l, [](double x) {
		return x;
	});
}

std::vector<Feature*> GLM::getFeatureList() {
	std::cerr << "GLM error: This operation is currently unsupported.";
	std::cerr << std::endl;
	throw std::exception();
}
