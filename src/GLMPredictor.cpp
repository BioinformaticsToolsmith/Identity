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
 * GLMPredictor.cpp
 *
 *  Created on: Jul 9, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#include "GLMPredictor.h"

GLMPredictor::GLMPredictor() {
	canDelete = false;
}

/**
 * If mode is true, calculate identity returns 0.0 or 1.0
 * Otherwise, it returns a number between 0.0 and 1.0
 */
GLMPredictor::GLMPredictor(std::vector<Feature*> &featList, bool mode) {
	if (featList.empty()) {
		std::cerr << "GLMPredictor error: Empty feature list. ";
		std::cerr << std::endl;
		throw std::exception();
	}

	isClassification = mode;

	// Get the bias and remove feature
	bias = featList.at(0)->getW();
	delete featList[0];
	featList.erase(featList.begin());
	featNum = featList.size();

	for (int u = 0; u < featList.size(); u++) {
		featList.at(u)->setTableIndex(u);
	}

	minList = new double[featNum];
	maxMinList = new double[featNum];
	distNum = 0;
	selectNum = 0;
	singleFeatNum = 0;
	for (int i = 0; i < featNum; i++) {
		auto f = featList.at(i);

		double min = f->getNormP1();
		double max = f->getNormP2();
		minList[i] = min;
		maxMinList[i] = max - min;

		if (f->getNumOfComp() == 0 && f->getName().compare("constant") != 0) {
			singleFeatNum++;
		}

		if (f->getIsDistance()) {
			distNum++;
		}

		if (f->getIsSelected()) {
			selectNum++;
		}
	}

	distIndexList = new int[distNum];
	wList = new double[selectNum];
	selectedIndexList = new int[selectNum];
	int x = 0;
	int w = 0;
	for (int i = 0; i < featNum; i++) {
		auto f = featList.at(i);

		if (f->getIsDistance()) {
			distIndexList[x] = i;
			x++;
		}

		if (f->getIsSelected()) {
			wList[w] = f->getW();
			selectedIndexList[w] = i;
			w++;
		}
	}

	// singleFeatNum = funIndexList.size();
	// funIndexArray = funIndexList.data();

	// Fill expansion list
	// The first part of this array is not filled
	expList = new std::pair<int, int>[featNum];
	for (int i = singleFeatNum; i < featNum; i++) {
		auto f = featList.at(i);

		if (f->getNumOfComp() == 1) {
			int c1 = f->getCompOneIndex();
			expList[i] = std::make_pair(c1, c1);
		}

		if (f->getNumOfComp() == 2) {
			int c1 = f->getCompOneIndex();
			int c2 = f->getCompTwoIndex();
			expList[i] = std::make_pair(c1, c2);
		}
	}

	canDelete = true;
}

GLMPredictor::GLMPredictor(const GLMPredictor &o) {
	copy(o);
}

GLMPredictor& GLMPredictor::operator=(const GLMPredictor &o) {
	copy(o);
	return *this;
}

void GLMPredictor::copy(const GLMPredictor &o) {
	canDelete = o.canDelete;
	isClassification = o.isClassification;
	singleFeatNum = o.singleFeatNum;
	featNum = o.featNum;
	bias = o.bias;
	distNum = o.distNum;
	selectNum = o.selectNum;

	minList = new double[featNum];
	maxMinList = new double[featNum];
	expList = new std::pair<int, int>[featNum];

	for (int i = 0; i < featNum; i++) {
		minList[i] = o.minList[i];
		maxMinList[i] = o.maxMinList[i];
		expList[i] = o.expList[i];
	}

	distIndexList = new int[distNum];
	for (int i = 0; i < distNum; i++) {
		distIndexList[i] = o.distIndexList[i];
	}

	wList = new double[selectNum];
	selectedIndexList = new int[selectNum];
	for (int i = 0; i < selectNum; i++) {
		wList[i] = o.wList[i];
		selectedIndexList[i] = o.selectedIndexList[i];
	}
}

GLMPredictor::~GLMPredictor() {
	if (canDelete) {
		delete[] minList;
		delete[] maxMinList;
		delete[] distIndexList;
		delete[] wList;
		delete[] expList;
		delete[] selectedIndexList;
	}
}

int GLMPredictor::getFeatNum() const {
	return featNum;
}

int GLMPredictor::getSingleFeatNum() const {
	return singleFeatNum;
}
