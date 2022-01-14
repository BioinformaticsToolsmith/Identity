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
 * IdentityCalculator.cpp
 *
 *     Created on: Dec 24, 2020
 *  Refactored on: Mar 13, 2021
 *         Author: Hani Zakaria Girgis, PhD
 *        Purpose: This class is all is needed to calculated identity scores
 *           Note: Each operation has its own thread number
 */

template<class V>
IdentityCalculator<V>::IdentityCalculator(DataGenerator *g, int threadNum,
		double t, bool skip, bool relax, std::string modelFile) {
	//g = generator;
	threshold = t;
	canSkip = skip;
	canRelax = relax;

	/**
	 * These data are needed for statisticians
	 */
	monoHistSize = Parameters::getAlphabetSize();

	double *temp = g->getCompositionList();
	compositionList = new double[monoHistSize];
	for (int i = 0; i < monoHistSize; i++) {
		compositionList[i] = temp[i];
	}

	k = g->getK();
	kHistSize = g->getHistogramSize();
	keyList = Util::makeKeyList(kHistSize, k);

	/**
	 * Train and prepare the predictor
	 */
	GLMRegressor regressor(g->getFeatures(), g->getLabels(), 0.0, threadNum,
			g->getK());
	regressor.start();

	absError = regressor.getAbsError();
	if (canRelax) {
		threshold -= regressor.getAbsError();
	}

	// Free memory used by the training and the validation data
	g->clearData();

	std::vector<Feature*> featList = regressor.getFeatureList(); // @suppress("Invalid arguments")
	featNum = featList.size() - 1; // The bias has not been removed yet.
	for (auto f : featList) {
		if (f->getNumOfComp() == 0 && f->getName().compare("constant") != 0) {
			funIndexList.push_back(f->getFunIndex());
		}
	}
	singleFeatNum = funIndexList.size();
	funIndexArray = funIndexList.data();

	// Write feature list
	if (!modelFile.empty()) {
		Serializer serializer(&featList, compositionList, k, kHistSize,
				absError, g->getMaxLength(), modelFile);
	}

	// This predictor removes the bias from the feature list
	p = GLMPredictor(featList, false);

	// featList is not needed beyond this point
	for (auto f : featList) {
		delete f;
	}
	featList.clear();

	// Construct table builders
	kTable = new KmerHistogram<uint64_t, V>(k);
	monoTable = new KmerHistogram<uint64_t, uint64_t>(1);
}

template<class V>
IdentityCalculator<V>::IdentityCalculator(Serializer &serializer, double t,
		bool skip, bool relax) {

	// Parse model
	// Serializer serializer(modelFile);

	threshold = t;
	canSkip = skip;
	canRelax = relax;

	/**
	 * These data are needed for statisticians
	 */
	monoHistSize = Parameters::getAlphabetSize();
	double *temp = serializer.getCompList();
	compositionList = new double[monoHistSize];
	for (int i = 0; i < monoHistSize; i++) {
		compositionList[i] = temp[i];
	}

	k = serializer.getK();
	kHistSize = serializer.getHistSize();
	absError = serializer.getAbsError();
	if (canRelax) {
		threshold -= absError;
	}

	keyList = Util::makeKeyList(kHistSize, k);

	vector<Feature*> *featList = serializer.getFeatList();
	featNum = featList->size() - 1; // The bias has not been removed yet.
	for (auto f : *featList) {
		if (f->getNumOfComp() == 0 && f->getName().compare("constant") != 0) {
			funIndexList.push_back(f->getFunIndex());
		}
	}
	singleFeatNum = funIndexList.size();
	funIndexArray = funIndexList.data();

	// This predictor removes the bias from the feature list
	p = GLMPredictor(*featList, false);

	// Construct table builders
	kTable = new KmerHistogram<uint64_t, V>(k);
	monoTable = new KmerHistogram<uint64_t, uint64_t>(1);
}

template<class V>
IdentityCalculator<V>::~IdentityCalculator() {
	delete[] compositionList;
	delete[] keyList;
	delete monoTable;
	delete kTable;
}

/**
 * One vs. all
 * Memory: The returned array is allocated on the heap. It is the
 * 	client responsibility to free its memory.
 */
template<class V>
double* IdentityCalculator<V>::score(V *kHist1, V **kHist2List,
		uint64_t *monoHist1, uint64_t **monoHist2List, int listSize,
		int threadNum, int len1, int *len2List) {

	double *v = new double[listSize];

#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < listSize; i++) {
		double ratio = calcRatio(len1, len2List[i]);
		if (canSkip && ratio < threshold) {
			v[i] = 0.0;
		} else {
			v[i] = score(kHist1, kHist2List[i], monoHist1, monoHist2List[i],
					ratio, len1, len2List[i]);
		}
	}
	return v;
}

/**
 * All vs. all in the same block filter too short and too long pairs
 */
template<class V>
Matrix IdentityCalculator<V>::score(V **kHistList, uint64_t **monoHistList,
		int listSize, int threadNum, int *lenList) {

	Matrix m(listSize, listSize, 0.0);

	for (int i = 0; i < listSize; i++) {
		m(i, i) = 1.0;
	}

	for (int i = 0; i < listSize; i++) {
#pragma omp parallel for schedule(static) num_threads(threadNum)
		for (int j = i + 1; j < listSize; j++) {
			double ratio = calcRatio(lenList[i], lenList[j]);
			if (!canSkip || ratio >= threshold) {
				double r = score(kHistList[i], kHistList[j], monoHistList[i],
						monoHistList[j], ratio, lenList[i], lenList[j]);
				m(i, j) = r;
				m(j, i) = r;
			}
		}
	}
	return m;
}

/**
 * Calculate all versus all
 * The result matrix is placed on the heap
 * It is the responsibility of the client to free allocated memory
 */
//template<class V>
//Matrix* IdentityCalculator<V>::scoreHeap(V **kHistList, uint64_t **monoHistList,
//		int listSize, int threadNum, int *lenList) {
//
//	Matrix *m = new Matrix(listSize, listSize, 1.0);
//	for (int i = 0; i < listSize; i++) {
//
//#pragma omp parallel for schedule(static) num_threads(threadNum)
//		for (int j = i + 1; j < listSize; j++) {
//			double ratio = calcRatio(lenList[i], lenList[j]);
//			if (canSkip && ratio < threshold) {
//				m->item(i, j) = 0.0;
//				m->item(j, i) = 0.0;
//			} else {
//				double r = score(kHistList[i], kHistList[j], monoHistList[i],
//						monoHistList[j], ratio);
//				m->item(i, j) = r;
//				m->item(j, i) = r;
//			}
//		}
//	}
//
//	return m;
//}
template<class V>
double IdentityCalculator<V>::getError() const {
	return absError;
}

template<class V>
int IdentityCalculator<V>::getK() const {
	return k;
}

/**
 * This method calculates the k-mer histograms and the mono
 * histograms. It frees memory used by the sequences
 * stored in the block.
 */
template<class V>
std::tuple<V**, uint64_t**, std::string**, int*> IdentityCalculator<V>::unpackBlock(
		Block *block, int threadNum) {
	int size = block->size();
	V **kHistList = new V*[size];
	uint64_t **monoHistList = new uint64_t*[size];
	std::string **infoList = new std::string*[size];
	int *lenList = new int[size];

#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < size; i++) {
		auto p = block->at(i);
		infoList[i] = p.first;
		std::string *seq = p.second;
		kHistList[i] = kTable->build(seq);
		monoHistList[i] = monoTable->build(seq);
		// A check
		if (Util::isAllZeros(kHistList[i], kHistSize)
				|| Util::isAllZeros(monoHistList[i], monoHistSize)) {
#pragma omp critical (unpackBlock)
			{
				std::cerr << "Found histograms made of all zeros: ";
				std::cerr << *infoList[i] << std::endl;
				throw std::exception();
			}
		}

		lenList[i] = seq->length();
		delete seq;
	}
	block->clear();
	delete block;

	return std::make_tuple(kHistList, monoHistList, infoList, lenList);
}

template<class V>
void IdentityCalculator<V>::freeBlock(
		std::tuple<V**, uint64_t**, std::string**, int*> t, int size,
		int threadNum) {
	auto kHistList = get < 0 > (t);
	auto monoHistList = get < 1 > (t);
	auto infoList = get < 2 > (t);
	auto lenList = get < 3 > (t);

#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < size; i++) {
		delete[] kHistList[i];
		delete[] monoHistList[i];
		delete infoList[i];
	}

	delete[] kHistList;
	delete[] monoHistList;
	delete[] infoList;
	delete[] lenList;
}

template<class V>
int IdentityCalculator<V>::getKHistSize() const {
	return kHistSize;
}

template<class V>
int IdentityCalculator<V>::getMonoHistSize() const {
	return monoHistSize;
}

template<class V>
bool IdentityCalculator<V>::getCanSkip() const {
	return canSkip;
}

template<class V>
void IdentityCalculator<V>::setCanSkip(bool yesNo) {
	canSkip = yesNo;
}
