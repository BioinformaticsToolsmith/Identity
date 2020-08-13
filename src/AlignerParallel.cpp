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
 * AlignerParallel.cpp
 *
 *  Created on: Jul 20, 2020
 *      Author: Hani Zakaria Girgis, PhD
 *  How to use it: (1) Initialize, (2) Set block A, and (3) Process block B.
 */

template<class V>
AlignerParallel<V>::AlignerParallel(int kmer, int hSize, double t,
		double e, // @suppress("Class members should be properly initialized")
		bool f, double *compList, ITransformer *transformer, std::string d,
		int tNum, std::string oFile) {
	k = kmer;
	histSize = hSize;
	threshold = t;
	// error = e;
	relaxThreshold = t - e;
	isLengthFilter = f;
	compositionList = compList;
	dlm = d;
	threadNum = tNum;

	out = std::ofstream(oFile.c_str(), std::ios::out);

	kTable = new KmerHistogram<uint64_t, V>(k);
	monoTable = new KmerHistogram<uint64_t, uint64_t>(1);

	keyList = new uint8_t[histSize * k];
	int alphaSize = Parameters::getAlphabetSize();
	for (int c = k - 1; c >= 0; c--) {
		for (int r = 0; r < histSize; r++) {
			keyList[(r * k) + k - 1 - c] = (r / ((uint64_t) pow(alphaSize, c)))
					% alphaSize;
		}
	}

	auto featList = transformer->getFeatureList();
	featNum = featList.size() - 1; // The bias has not been removed yet.
	for (auto f : featList) {
		if (f->getNumOfComp() == 0 && f->getName().compare("constant") != 0) {
			funIndexList.push_back(f->getFunIndex());
		}
	}
	singleFeatNum = funIndexList.size();
	funIndexArray = funIndexList.data();

	// This predictor removes the bias from the feature list
	predictor = GLMPredictor(featList, false);

	// featList is not needed beyond this point
	for (auto f : featList) {
		delete f;
	}
	featList.clear();

	isInitialized = false;
}

template<class V>
AlignerParallel<V>::~AlignerParallel() {
	out.close();
	if (isInitialized) {
		clearAMemory(kHistList, monoHistList, infoList, lenList, sizeA);
	}

	delete kTable;
	delete monoTable;
	delete[] keyList;
}

/**
 * This method clears memory used by block A
 */
template<class V>
void AlignerParallel<V>::clearAMemory(V **kHistList, uint64_t **monoHistList,
		std::string **infoList, int *lenList, int size) {
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

/**
 * Note! This method deallocates almost all memory
 * held by the block---except it keeps sequence informations.
 */
template<class V>
void AlignerParallel<V>::setBlockA(Block *block, bool isAllVsAll) {
	if (isInitialized) {
		clearAMemory(kHistList, monoHistList, infoList, lenList, sizeA);
	} else {
		isInitialized = true;
	}

	sizeA = block->size();
	auto tup = unpackBlock(block);
	kHistList = std::get < 0 > (tup);
	monoHistList = std::get < 1 > (tup);
	infoList = std::get < 2 > (tup);
	lenList = std::get < 3 > (tup);

	if (isAllVsAll) {
		std::future<void> printTask;
		for (int i = 0; i < sizeA; i++) {
			auto printList = makeEmptyResult(threadNum,
					((sizeA - i - 1) / threadNum) + 1);

#pragma omp parallel for schedule(static) num_threads(threadNum)
			for (int j = i + 1; j < sizeA - 1; j++) {
				double minimum = lenList[i];
				double maximum = lenList[j];
				if (maximum < minimum) {
					minimum = lenList[j];
					maximum = lenList[i];
				}

				if (isLengthFilter && (minimum / maximum < threshold)) {
					continue;
				}

				Statistician < V
						> s(histSize, k, kHistList[i], kHistList[j],
								monoHistList[i], monoHistList[j],
								compositionList, keyList);
				double data[featNum];
				s.calculate(funIndexArray, singleFeatNum, data);
				double res = predictor.calculateIdentity(data);

				if (res >= relaxThreshold) {
					printList->at(omp_get_thread_num())->push_back(
							std::make_pair(infoList[j], res));
				}
			}

			// Wait for the previous print task to finish
			if (printTask.valid()) {
				printTask.get();
				threadNum++;
			}

			// Do not run a writing task if there is nothing to be written
			for (auto ptr : *printList) {
				if (ptr->size() > 0) {
					threadNum--;
					printTask = std::async([this, printList, i]() {
						output(printList, infoList[i]);
					});
					break;
				}
			}
		}

		// Wait for the last print case
		if (printTask.valid()) {
			printTask.get();
			threadNum++;
		}
	}
}

/**
 * This method calculate the k-mer histograms and the mono
 * histograms. It frees memory used by the sequences
 * stored in the block.
 */
template<class V>
std::tuple<V**, uint64_t**, std::string**, int*> AlignerParallel<V>::unpackBlock(
		Block *block) {
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
		lenList[i] = seq->length();
		delete seq;
	}
	delete block;

	return std::make_tuple(kHistList, monoHistList, infoList, lenList);
}

template<class V>
Result* AlignerParallel<V>::makeEmptyResult(int outerSize, int innerSize) {
	auto *printList = new Result();
	printList->reserve(outerSize);
	for (int e = 0; e < outerSize; e++) {
		auto v = new std::vector<pair<std::string*, double> >();
		v->reserve(innerSize);
		printList->push_back(v);
	}
	return printList;
}

/**
 * Print some results to a file
 * Memory: This method frees the memory used by the print list.
 */
template<class V>
void AlignerParallel<V>::output(
		std::vector<std::vector<pair<std::string*, double> >*> *printList,
		string *info1) {
	for (auto v : *printList) {
		for (auto p : *v) {
			double res = p.second;

			if (res > 1.0) {
				res = 1.0;
			} else if (res < 0.0) {
				res = 0.0;
			}

			out << *info1 << dlm << *p.first << dlm << std::setprecision(4)
					<< res << std::endl;
		}
		delete v;
	}
	delete printList;
}

/**
 * Note! Memory allocated to the block are freed here.
 */
template<class V>
void AlignerParallel<V>::processBlockB(Block *block) {
	int sizeB = block->size();
	auto tup = unpackBlock(block);
	auto kHistListB = std::get < 0 > (tup);
	auto monoHistListB = std::get < 1 > (tup);
	auto infoListB = std::get < 2 > (tup);
	auto lenListB = std::get < 3 > (tup);

	std::future<void> printTask;
	for (int i = 0; i < sizeA; i++) {
		auto printList = makeEmptyResult(threadNum, (sizeB / threadNum) + 1);

#pragma omp parallel for schedule(static) num_threads(threadNum)
		for (int h = 0; h < sizeB; h++) {

			double minimum = lenList[i];
			double maximum = lenListB[h];
			if (maximum < minimum) {
				minimum = lenListB[h];
				maximum = lenList[i];
			}

			if (isLengthFilter && (minimum / maximum < threshold)) {
				continue;
			}

			Statistician < V
					> s(histSize, k, kHistList[i], kHistListB[h],
							monoHistList[i], monoHistListB[h], compositionList,
							keyList);
			double data[featNum];
			s.calculate(funIndexArray, singleFeatNum, data);
			double res = predictor.calculateIdentity(data);

			if (res >= relaxThreshold) {
				printList->at(omp_get_thread_num())->push_back(
						std::make_pair(infoListB[h], res));
			}
		}

		// Wait for the previous print task to finish
		if (printTask.valid()) {
			printTask.get();
			threadNum++;
		}

		// Do not run a writing task if there is nothing to be written
		for (auto ptr : *printList) {
			if (ptr->size() > 0) {
				threadNum--;
				printTask = std::async([this, printList, i]() {
					output(printList, infoList[i]);
				});
				break;
			}
		}
	}

	// Wait for the last print case
	if (printTask.valid()) {
		printTask.get();
		threadNum++;
	}

	clearAMemory(kHistListB, monoHistListB, infoListB, lenListB, sizeB);
}

template<class V>
int AlignerParallel<V>::getThreadNum() const {
	return threadNum;
}

template<class V>
void AlignerParallel<V>::setThreadNum(int threadNum) {
	this->threadNum = threadNum;
}
