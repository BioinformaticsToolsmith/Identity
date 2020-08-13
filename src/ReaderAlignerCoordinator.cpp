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
 * ReaderAlignerCoordinator.cpp
 *
 *  Created on: Nov 15, 2019
 *      Author: Hani Z. Girgis, PhD
 */

#include "ReaderAlignerCoordinator.h"

ReaderAlignerCoordinator::ReaderAlignerCoordinator(int workerNumIn,
		int blockSizeIn, char m, double t, bool r) {
	workerNum = workerNumIn;
	blockSize = blockSizeIn;
	mode = m;
	threshold = t;
	canRelax = r;
}

ReaderAlignerCoordinator::~ReaderAlignerCoordinator() {

}

void ReaderAlignerCoordinator::alignAllVsAll(string fileIn, string fileOut,
		string dlm) {
	alignFileVsFile1(fileIn, fileIn, fileOut, dlm, true);
}

void ReaderAlignerCoordinator::alignQueryVsAll(string fileDb, string fileQry,
		string fileOut, string dlm) {
	alignFileVsFile2(fileDb, fileQry, fileOut, dlm, false);
}

void ReaderAlignerCoordinator::alignFileVsFile1(string fileDb, string fileQry,
		string fileOut, string dlm, bool isAllVsAll) {

	DataGenerator *g;
	if (isAllVsAll) {
		g = new SynDataGenerator(fileDb, threshold, workerNum);
	} else {
		g = new SynDataGenerator(fileDb, fileQry, threshold, workerNum);
	}

	int64_t maxLength = g->getMaxLength();

	// Determine histogram data type
	if (maxLength <= std::numeric_limits<int8_t>::max()) {
		helper<int8_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g);
	} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
		helper<int16_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g);
	} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
		helper<int32_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g);
	} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
		helper<int64_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g);
	} else {
		std::cout << "ReaderAlignerCoordinator warning: ";
		std::cout << "Overflow is possible however unlikely.";
		std::cout << std::endl;
		std::cout << "A histogram entry is 64 bits." << std::endl;
		helper<int64_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g);
	}
}

template<class V>
void ReaderAlignerCoordinator::helper(string fileDb, string fileQry,
		string fileOut, string dlm, bool isAllVsAll, DataGenerator *g) {

	int hSize = g->getHistogramSize();
	int k = g->getK();

	ITransformer *t;
	double error = 0.0;
	switch (mode) {
	case C: {
		std::cout << "Mode is classification." << std::endl;
		GLMClassifier *c = new GLMClassifier(g->getFeatures(), g->getLabels(),
				threshold, workerNum, k);
		c->start();
		t = c;
	}
		break;
	case R: {
		std::cout << "Mode is regression." << std::endl;
		// If regression alone, it must learn the entire function not part of it.
		// That is why no threshold, i.e. threshold = 0.
		GLMRegressor *r = new GLMRegressor(g->getFeatures(), g->getLabels(),
				0.0, workerNum, k);
		r->start();
		error = r->getAbsError();
		t = r;
	}
		break;
	default:
		std::cerr << "ReaderAlignerCoordinator error: Invalid mode.";
		std::cerr << std::endl;
		throw std::exception();
	}

	if (canRelax) {
		std::cout << "Relaxing the threshold by " << error << std::endl;
	} else {
		error = 0.0;
	}

	// Free memory used by the training and the validation data
	g->clearData();

	std::cout
			<< "Calculating the identity scores. This step may take long time ..."
			<< std::endl;

	FastaReader qryReader(fileQry, blockSize);
	AlignerParallel<V> aligner(k, hSize, threshold, error, true,
			g->getCompositionList(), t, dlm, workerNum, fileOut);
	if (isAllVsAll) {
		// Process the first block versus itself.
		aligner.setBlockA(qryReader.read(), isAllVsAll);
		while (qryReader.isStillReading()) {
			// Construct a database reader
			FastaReader dbReader(fileDb, blockSize, qryReader.getCurrentPos(),
					qryReader.getMaxLen());

			// Make sure there is one free thread for reading
			aligner.setThreadNum(workerNum - 1);
			// Start a reading task
			LockFreeQueue<Block*, 1000> buffer;
			auto readFuture = std::async([&dbReader, &buffer]() -> int {
				int blockRead = 0;
				while (dbReader.isStillReading()) {
					buffer.push(dbReader.read());
					blockRead++;
				}
				return blockRead;
			});

			// Align each query block versus this database block
			int blockProcessed = 0;
			int blockRead = -1;
			while (true) {
				if (buffer.size() > 0) {
					if (!dbReader.isStillReading()) {
						// The reading thread is done. Use it in the aligner.
						aligner.setThreadNum(workerNum);
					}

					aligner.processBlockB(buffer.front());
					buffer.pop();
					blockProcessed++;

				} else if (!dbReader.isStillReading()) {
					// The reading thread is done. Use it in the aligner.
					aligner.setThreadNum(workerNum);

					if (readFuture.valid()) {
						blockRead = readFuture.get();
					}
					if (blockProcessed == blockRead) {
						break;
					}
				}
			}

			aligner.setBlockA(qryReader.read(), isAllVsAll);
		}
	} else {
		while (qryReader.isStillReading()) {
			aligner.setBlockA(qryReader.read(), isAllVsAll);

			// Construct a database reader
			FastaReader dbReader(fileDb, blockSize);

			// Make sure there is one free thread for reading
			aligner.setThreadNum(workerNum - 1);

			// Start a reading task
			LockFreeQueue<Block*, 1000> buffer;
			auto readFuture = std::async([&dbReader, &buffer]() -> int {
				int blockRead = 0;
				while (dbReader.isStillReading()) {
					buffer.push(dbReader.read());
					blockRead++;
				}
				return blockRead;
			});

			// Align each query block versus this database block
			int blockProcessed = 0;
			int blockRead = -1;
			while (true) {
				if (buffer.size() > 0) {
					if (!dbReader.isStillReading()) {
						// The reading thread is done. Use it in the aligner.
						aligner.setThreadNum(workerNum);
					}

					aligner.processBlockB(buffer.front());
					buffer.pop();

					blockProcessed++;
				} else if (!dbReader.isStillReading()) {
					// The reading thread is done. Use it in the aligner.
					aligner.setThreadNum(workerNum);

					if (readFuture.valid()) {
						blockRead = readFuture.get();
					}

					if (blockProcessed == blockRead) {
						break;
					}
				}
			}
		}
	}

	delete g;
	delete t;
}

void ReaderAlignerCoordinator::alignFileVsFile2(string fileDb, string fileQry,
		string fileOut, string dlm, bool isAllVsAll) {
	DataGenerator *g;
	if (isAllVsAll) {
		g = new SynDataGenerator(fileDb, threshold, workerNum);
	} else {
		g = new SynDataGenerator(fileDb, fileQry, threshold, workerNum);
	}

	// Make the keyList in digital format once
	int alphaSize = Parameters::getAlphabetSize();
	int histogramSize = g->getHistogramSize();
	int k = g->getK();
	uint8_t keyList[histogramSize * k];
	for (int c = k - 1; c >= 0; c--) {
		for (int r = 0; r < histogramSize; r++) {
			keyList[(r * k) + k - 1 - c] = (r / ((uint64_t) pow(alphaSize, c)))
					% alphaSize;
		}
	}

	ITransformer *t;
	double error = 0.0;
	switch (mode) {
	case C: {
		std::cout << "Mode is classification." << std::endl;
		GLMClassifier *c = new GLMClassifier(g->getFeatures(), g->getLabels(),
				threshold, workerNum, k);
		c->start();
		t = c;
	}
		break;
	case R: {
		std::cout << "Mode is regression." << std::endl;
		// If regression alone, it must learn the entire function not part of it.
		// That is why no threshold, i.e. threshold = 0.
		GLMRegressor *r = new GLMRegressor(g->getFeatures(), g->getLabels(),
				0.0, workerNum, k);
		r->start();
		error = r->getAbsError();
		t = r;
	}
		break;
	default:
		std::cerr << "ReaderAlignerCoordinator error: Invalid mode.";
		std::cerr << std::endl;
		throw std::exception();
	}

	if (canRelax) {
		std::cout << "Relaxing the threshold by " << error << std::endl;
	} else {
		error = 0.0;
	}

	// Free memory used by the training and the validation data
	g->clearData();

	// Make sure you have one thread for reading
	workerNum--;

	std::cout
			<< "Calculating the identity scores. This step may take long time ..."
			<< std::endl;

	if (!isAllVsAll) {
		string temp(fileDb);
		fileDb = fileQry;
		fileQry = temp;
	}

	// Construct a database reader.
	FastaReader dbReader(fileDb, blockSize);

	// Open output file
	std::ofstream out(fileOut.c_str(), std::ios::out);

	while (dbReader.isStillReading()) {
		// Construct a query reader
		FastaReader *qryReader;
		if (isAllVsAll) {
			qryReader = new FastaReader(fileQry, blockSize,
					dbReader.getCurrentPos(), dbReader.getMaxLen());
		} else {
			qryReader = new FastaReader(fileQry, blockSize);
		}

		// Start a reading task
		LockFreeQueue<pair<Block*, bool>, 1000> buffer;
		auto readFuture = std::async([qryReader, isAllVsAll, &buffer]() {
			// When this boolean is true, the db first block and the qry first block
			// are the same.
			bool isQryFirst = isAllVsAll;
			while (qryReader->isStillReading()) {
				// Read a query block
				Block *qryBlock = qryReader->read();
				buffer.push(make_pair(qryBlock, isQryFirst));

				if (isQryFirst) {
					// ToDo: Review this regarding one versus all
					isQryFirst = false;
				}
			}
		});

		// Read a database block.
		Block *dbBlock = dbReader.read();

		// Start concurrent aligner tasks
		vector<Aligner*> alignerList;
		alignerList.reserve(workerNum);
		vector<future<std::pair<bool, stringstream*> > > futureList;
		futureList.reserve(workerNum);
		for (int i = 0; i < workerNum; i++) {
			Aligner *aligner = new Aligner(g, t, dbBlock, dlm, true, threshold,
					error, keyList);
			alignerList.push_back(aligner);
			futureList.push_back(
					std::async([aligner]() -> std::pair<bool, stringstream*> {
						return aligner->start();
					}));
		}

		// Read a block and pass it to one of the workers
		int nextIndex = 0;
		while (true) {
			if (buffer.size() > 0) {
				alignerList.at(nextIndex)->enqueueBlock(buffer.front());
				buffer.pop();
				nextIndex = (nextIndex + 1) % workerNum;
			} else if (!qryReader->isStillReading()) {
				break;
			} else {
				std::this_thread::sleep_for(std::chrono::seconds(1));
			}
		}

		// Tell the workers that no more blocks will be passed to them.
		for (int i = 0; i < workerNum; i++) {
			alignerList.at(i)->stop();
		}

		// Wait until all of the workers are done.
		for (int i = 0; i < workerNum; i++) {
			auto p = futureList.at(i).get();
			if (p.first) {
				out << p.second->rdbuf();
			}
		}

		// Free resources
		FastaReader::deleteBlock(dbBlock);

		// Delete the aligners
		for (int i = 0; i < workerNum; i++) {
			delete alignerList.at(i);
		}
		alignerList.clear();
		futureList.clear();
		delete qryReader;
	}
	cout << endl;

	// Close output file.
	out.close();

	delete g;
	delete t;
}
