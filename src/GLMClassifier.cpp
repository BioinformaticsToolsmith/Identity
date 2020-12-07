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
 * GLMClassifier.cpp
 *
 *  Created on: Mar 31, 2020
 *      Author: Dr. Hani Z. Girgis
 *
 */
#include "GLMClassifier.h"

GLMClassifier::GLMClassifier(const Matrix *f, const Matrix *l, double t, int c,
		int m, double b) :
		fModelTable(f), lModelTable(l), threshold(t), threadNum(c), minFeat(m), balance(
				b) {
	pipe = new std::vector<ITransformer*>();

	auto info = StatisticInfo::getInstance();
	std::vector<Feature*> *oList = info->getList();
	// Make a copy of the feature list
	fList = new std::vector<Feature*>();
	int s = oList->size();
	fList->reserve(s);
	for (int i = 0; i < s; i++) {
		auto f = new Feature(*oList->at(i));
		f->setTableIndex(i);
		fList->push_back(f);
	}
}

GLMClassifier::~GLMClassifier() {
	// Training tables are deleted right after training
	// Validation tables are deleted right after validation

	for (auto ptr : *pipe) {
		delete ptr;
	}
	pipe->clear();
	delete pipe;

	for (auto ptr : *fList) {
		delete ptr;
	}
	fList->clear();
	delete fList;

	clean(f5);
}

void GLMClassifier::start() {
	prepareData();
	train();
	validate();
}

/**
 * Fill the four matrices here
 * Balanced data: Positive examples = negative examples
 */
void GLMClassifier::prepareData() {
	std::cout << "Preparing data ..." << std::endl;
	// Determine matrices sizes
	int numPstv = 0;
	int numNgtv = 0;
	int numRow = lModelTable->getNumRow();

	for (int i = 0; i < numRow; i++) {
		if (lModelTable->item(i, 0) >= threshold) {
			numPstv++;
		} else {
			numNgtv++;
		}
	}

	// Make sure we have positive and negative examples
	if (numPstv == 0) {
		std::cerr << "GLMClassifier error: No positives." << std::endl;
		throw std::exception();
	}
	if (numNgtv == 0) {
		std::cerr << "GLMClassifier error: No negatives." << std::endl;
		throw std::exception();
	}

	// I need numPstv and numNgtv to be even
	if (numPstv % 2 != 0) {
		numPstv--;
	}
	if (numNgtv % 2 != 0) {
		numNgtv--;
	}

	// Build two balanced sets for training and validation
	int p = 0; // Number of collected training and validation positive examples
	int n = 0; // Number of collected training and validation negative examples
	int t = 0; // Number of collected training and validation examples
	int v = 0; // Number of collected training and validation examples
	int i = 0; // Index in model table

	int s = std::min(numPstv, numNgtv);
	int pstvS = s;
	int ngtvS = balance * s;
	if (ngtvS > numNgtv) {
		ngtvS = numNgtv;
	}

	int tableS = (pstvS + ngtvS) / 2;
	int trainList[tableS]; // Array holding training indexes
	int validateList[tableS]; // Array holding validation indexes

	// Label matrices
	lTrainTable = new Matrix(tableS, 1, 1.0);
	lValidateTable = new Matrix(tableS, 1, 1.0);

	while (p < pstvS || n < ngtvS) {
		double id = lModelTable->at(i, 0);
		if (id >= threshold && p < pstvS) {
			if (p % 2 == 0) {
				//lTrainTable->at(t, 0) = id;
				trainList[t] = i;
				t++;
			} else {
				//lValidateTable->at(v, 0) = id;
				validateList[v] = i;
				v++;
			}
			p++;
		} else if (id < threshold && n < ngtvS) {
			if (n % 2 == 0) {
				lTrainTable->at(t, 0) = 0.0;
				//lTrainTable->at(t, 0) = id;
				trainList[t] = i;
				t++;
			} else {
				lValidateTable->at(v, 0) = 0.0;
				// lValidateTable->at(v, 0) = id;
				validateList[v] = i;
				v++;
			}
			n++;
		}
		i++;
	}

	// Feature matrices
	fTrainTable = fModelTable->subMatrix(trainList, tableS);
	fValidateTable = fModelTable->subMatrix(validateList, tableS);

	std::cout << "\tSimilar pair count: " << p << std::endl;
	std::cout << "\tDissimilar pair count: " << n << std::endl;
	std::cout << "\tTraining size: " << t << std::endl;
	std::cout << "\tValidation size: " << v << std::endl;
}

pair<Matrix, std::vector<Feature*> > GLMClassifier::selectFeatures(Matrix &t4,
		std::vector<Feature*> &f4) {
	// Better if results in 0.1% increase in accuracy
	auto isNewBetter = [](double newV, double oldV) {
		return (newV - oldV > 0.001) ? true : false;
	};
	BestFirst<GLM> selector(t4, *lTrainTable, f4, GLM::classifierFactory,
			Evaluator::acc, isNewBetter, true, threadNum, minFeat);
	return std::make_pair(selector.transform(t4), selector.getFeatureList());
}

pair<Matrix, GLM*> GLMClassifier::trainGLM(Matrix &t5) {
	GLM *glm = GLM::classifierFactoryHeap(t5, *lTrainTable);
	Matrix t6 = glm->transform(t5);
	return std::make_pair(t6, glm);
}

void GLMClassifier::evaluate(Matrix &o, Matrix &p) {
	/**
	 * Those metrics could be due to training or validation
	 * They are due to validation if the validation step is done
	 */
	acc = Evaluator::acc(o, p);
	sens = Evaluator::sens(o, p);
	spec = Evaluator::spec(o, p);
	std::cout << "\tAccuracy: " << acc << std::endl;
	std::cout << "\tSensitivity: " << sens << std::endl;
	std::cout << "\tSpecificity: " << spec << std::endl;
}

void GLMClassifier::clean(std::vector<Feature*> &list) {
	for (auto f : list) {
		delete f;
	}
}

/**
 * Train and optimize the pipeline
 */
void GLMClassifier::train() {
	Normalizer normalizer1(*fTrainTable, *fList);
	Matrix t1 = normalizer1.transform(*fTrainTable);
	std::vector<Feature*> f1 = normalizer1.getFeatureList();

	SimConverter sim(t1, f1);
	Matrix t2 = sim.transform(t1);
	std::vector<Feature*> f2 = sim.getFeatureList();
	clean(f1);

	FeatureExpander expander(t2, f2);
	Matrix t3 = expander.transform(t2);
	std::vector<Feature*> f3 = expander.getFeatureList();
	clean(f2);

	Normalizer normalizer2(t3, f3);
	Matrix t4 = normalizer2.transform(t3);
	std::vector<Feature*> f4 = normalizer2.getFeatureList();
	clean(f3);

	auto p5 = selectFeatures(t4, f4);

	// auto p5 = selectFeatures(t2, f2); // No squared or paired statistics

	auto t5 = p5.first;
	f5 = p5.second;
	// clean(f2); // No squared or paired statistics
	clean(f4);

	auto p6 = trainGLM(t5);
	auto t6 = p6.first;
	auto glm = p6.second;

	std::cout << "Finished training." << std::endl;

	evaluate(*lTrainTable, t6);

	std::cout << "Optimizing ..." << std::endl;

	// List of single features
	std::vector<Feature*> singleList;
	// List of single indexes in the original table that includes all features
	int i = 0;
	for (auto ptr : f5) {
		if (ptr->getNumOfComp() == 0) {
			indexList.push_back(ptr->getTableIndex());

			auto f = new Feature(*ptr);
			f->setIsNormalized(false);
			f->setIsConverted(false);
			f->setTableIndex(i);
			singleList.push_back(f);
			i++;
		}
	}

	pipe->push_back(new Normalizer(singleList));
	pipe->push_back(new SimConverter(singleList));
	clean(singleList);

	for (int h = 0; h < f5.size(); h++) {
		auto ptr = f5[h];
		if (ptr->getNumOfComp() != 0) {
			ptr->setIsNormalized(false);
		}
		ptr->setTableIndex(h);
	}

	pipe->push_back(new FeatureExpander(f5));
	pipe->push_back(new Normalizer(f5));

	for (int h = 0; h < f5.size(); h++) {
		auto ptr = f5[h];
		if (ptr->getNumOfComp() != 0) {
			ptr->setIsNormalized(true);
		}
	}

	pipe->push_back(new BestFirst<GLM>(f5));
	pipe->push_back(glm);

	// For now. If it works move to GLM
	// Add weights to the features
	// The bias is the very first feature and it is not selected
	Feature *bias = new Feature(-1, "constant", true);
	Matrix w = glm->getWeights();

	bias->setW(w.at(0, 0));
	f5.insert(f5.begin(), bias);
	int g = 1;
	for (auto f : f5) {
		if (f->getIsSelected()) {
			f->setW(w.at(g, 0));
			g++;
		}
	}

	for (int u = 0; u < f5.size(); u++) {
		f5.at(u)->setTableIndex(u);
	}
	// End adding weights to the features

	// Clean up
	delete fTrainTable;
	delete lTrainTable;
}

void GLMClassifier::validate() {
	std::cout << "Validating ..." << std::endl;
	Matrix sub = fValidateTable->subMatrixByCol(indexList.data(),
			indexList.size());
	Matrix p = transform(sub);
	evaluate(*lValidateTable, p);

	// Clean up
	delete fValidateTable;
	delete lValidateTable;
}

/**
 * Runs the pipeline
 */
Matrix GLMClassifier::transform(const Matrix &f) {
	int s = pipe->size();

	if (s == 0) {
		std::cerr << "GLMClassifier error: Pipe is empty." << std::endl;
		throw std::exception();
	}

	Matrix m = pipe->at(0)->transform(f);
	for (int i = 1; i < s; i++) {
		m = pipe->at(i)->transform(m);
	}
	return m;
}

std::vector<Feature*> GLMClassifier::getFeatureList() {
	auto temp = copy(f5);
	for (auto f : temp) {
		f->setIsNormalized(false);
		f->setIsConverted(false);
	}

	return temp;
}

double GLMClassifier::getAcc() const {
	return acc;
}

double GLMClassifier::getSens() const {
	return sens;
}

double GLMClassifier::getSpec() const {
	return spec;
}
