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
 * ReaderAlignerCoordinator.cpp
 *
 *  Created on: Nov 15, 2019
 *      Author: Hani Z. Girgis, PhD
 */

#include "ReaderAlignerCoordinator.h"

ReaderAlignerCoordinator::ReaderAlignerCoordinator(
		int workerNumIn, // @suppress("Class members should be properly initialized")
		int blockSizeIn, double t, bool r, bool a, bool s, bool f,
		std::string file) {
	workerNum = workerNumIn;
	blockSize = blockSizeIn;
	threshold = t;
	canRelax = r;
	canReportAll = a;
	canSaveModel = s;
	canFillModel = f;
	modelFile = file;
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

	if (canFillModel) {
		Serializer serializer(modelFile);
		int64_t maxLength = serializer.getMaxLength();
		int hSize = serializer.getHistSize();
		int k = serializer.getK();
		double error = serializer.getAbsError();
		// Determine histogram data type
		if (maxLength <= std::numeric_limits<int8_t>::max()) {
			AlignerParallel<int8_t> aligner(serializer, threshold, canReportAll,
					dlm, workerNum, fileOut);
			helper1<int8_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, aligner);
		} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
			AlignerParallel<int16_t> aligner(serializer, threshold,
					canReportAll, dlm, workerNum, fileOut);
			helper1<int16_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll,
					aligner);
		} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
			AlignerParallel<int32_t> aligner(serializer, threshold,
					canReportAll, dlm, workerNum, fileOut);
			helper1<int32_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll,
					aligner);
		} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
			AlignerParallel<int64_t> aligner(serializer, threshold,
					canReportAll, dlm, workerNum, fileOut);
			helper1<int64_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll,
					aligner);
		} else {
			std::cout << "ReaderAlignerCoordinator warning: ";
			std::cout << "Overflow is possible however unlikely.";
			std::cout << std::endl;
			std::cout << "A histogram entry is 64 bits." << std::endl;

			AlignerParallel<int64_t> aligner(serializer, threshold,
					canReportAll, dlm, workerNum, fileOut);
			helper1<int64_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll,
					aligner);
		}
	} else {
		DataGenerator *g = nullptr;
		if (isAllVsAll) {
			g = new SynDataGenerator(fileDb, threshold, workerNum);
		} else {
			g = new SynDataGenerator(fileDb, fileQry, threshold, workerNum);
		}
		int64_t maxLength = g->getMaxLength();
		int hSize = g->getHistogramSize();
		int k = g->getK();
		// If regression alone, it must learn the entire function not part of it.
		// That is why no threshold, i.e. threshold = 0.
		GLMRegressor r(g->getFeatures(), g->getLabels(), 0.0, workerNum, k);
		r.start();
		double error = r.getAbsError();
		// Free memory used by the training and the validation data
		g->clearData();
		if (canRelax) {
			std::cout << "Relaxing the threshold by " << error << std::endl;
		} else {
			error = 0.0;
		}

		ITransformer *t = &r;
		// Determine histogram data type
		if (maxLength <= std::numeric_limits<int8_t>::max()) {
			AlignerParallel<int8_t> aligner(k, hSize, threshold, error,
					canReportAll, g->getCompositionList(), t, dlm, workerNum,
					maxLength, fileOut, modelFile);
			helper1<int8_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, aligner);
		} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
			AlignerParallel<int16_t> aligner(k, hSize, threshold, error,
					canReportAll, g->getCompositionList(), t, dlm, workerNum,
					maxLength, fileOut, modelFile);
			helper1<int16_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll,
					aligner);
		} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
			AlignerParallel<int32_t> aligner(k, hSize, threshold, error,
					canReportAll, g->getCompositionList(), t, dlm, workerNum,
					maxLength, fileOut, modelFile);
			helper1<int32_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll,
					aligner);
		} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
			AlignerParallel<int64_t> aligner(k, hSize, threshold, error,
					canReportAll, g->getCompositionList(), t, dlm, workerNum,
					maxLength, fileOut, modelFile);
			helper1<int64_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll,
					aligner);
		} else {
			std::cout << "ReaderAlignerCoordinator warning: ";
			std::cout << "Overflow is possible however unlikely.";
			std::cout << std::endl;
			std::cout << "A histogram entry is 64 bits." << std::endl;

			AlignerParallel<int64_t> aligner(k, hSize, threshold, error,
					canReportAll, g->getCompositionList(), t, dlm, workerNum,
					maxLength, fileOut, modelFile);
			helper1<int64_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll,
					aligner);
		}

		delete g;
	}
}

template<class V>
void ReaderAlignerCoordinator::helper1(string fileDb, string fileQry,
		string fileOut, string dlm, bool isAllVsAll,
		AlignerParallel<V> &aligner) {

	/**
	 * If the output file is not provided, identity is used in the training mode only
	 */
	if (fileOut.empty()) {
		std::cout << "Training is done." << std::endl;
		return;
	}

	std::cout
			<< "Calculating the identity scores. This step may take long time ..."
			<< std::endl;

	FastaReader qryReader(fileQry, blockSize);

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
}

void ReaderAlignerCoordinator::alignFileVsFile2(string fileDb, string fileQry,
		string fileOut, string dlm, bool isAllVsAll) {
	DataGenerator *g = nullptr;
	Serializer *serializer = nullptr;
	int64_t maxLength;
	if (canFillModel) {
		serializer = new Serializer(modelFile);
		maxLength = serializer->getMaxLength();
	} else {
		if (isAllVsAll) {
			g = new SynDataGenerator(fileDb, threshold, workerNum);
		} else {
			g = new SynDataGenerator(fileDb, fileQry, threshold, workerNum);
		}
		maxLength = g->getMaxLength();
	}

	// Determine histogram data type
	if (maxLength <= std::numeric_limits<int8_t>::max()) {
		helper2<int8_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g,
				serializer);
	} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
		helper2<int16_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g,
				serializer);
	} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
		helper2<int32_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g,
				serializer);
	} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
		helper2<int64_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g,
				serializer);
	} else {
		std::cout << "ReaderAlignerCoordinator warning: ";
		std::cout << "Overflow is possible however unlikely.";
		std::cout << std::endl;
		std::cout << "A histogram entry is 64 bits." << std::endl;
		helper2<int64_t>(fileDb, fileQry, fileOut, dlm, isAllVsAll, g,
				serializer);
	}

	if (g != nullptr) {
		delete g;
	}
	if (serializer != nullptr) {
		delete serializer;
	}
}

/**
 * Not for all versus all
 * Simple, wrote it to avoid the bug that showed when running identity on a small number of sequences
 */
template<class V>
void ReaderAlignerCoordinator::helper2_simple(string fileDb, string fileQry,
		string fileOut, string dlm, bool isAllVsAll, DataGenerator *g,
		Serializer *serializer) {

	bool canSkip = !canReportAll;

	IdentityCalculator<V> *id;
	if (canFillModel) {
		id = new IdentityCalculator<V>(*serializer, threshold, canSkip,
				canRelax);
	} else if (canSaveModel) {
		id = new IdentityCalculator<V>(g, workerNum, threshold, canSkip,
				canRelax, modelFile);
	} else {
		id = new IdentityCalculator<V>(g, workerNum, threshold, canSkip,
				canRelax);
	}

	/**
	 * If the output file is not provided, identity is used in the training mode only
	 */
	if (fileOut.empty()) {
		std::cout << "Training is done." << std::endl;
		return;
	}

	if (canRelax) {
		std::cout << "Relaxing the threshold" << std::endl;
	}

	std::cout << "Calculating the identity scores. ";
	std::cout << "This step may take long time ..." << std::endl;

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
		// Read a database block.
		Block *dbBlock = dbReader.read();
		int dbSize = dbBlock->size();
		auto dbTuble = id->unpackBlock(dbBlock, workerNum);
		V **dbKHistList = get<0>(dbTuble);
		uint64_t **dbMonoHistList = get<1>(dbTuble);
		std::string **dbInfoList = get<2>(dbTuble);
		int *dbLenList = get<3>(dbTuble);

		// Construct a query reader
		FastaReader qryReader(fileQry, blockSize);
		while (qryReader.isStillReading()) {
			// Read a query block
			Block *qryBlock = qryReader.read();		// Destroyed by the aligner
			int qrySize = qryBlock->size();
			auto qryTuble = id->unpackBlock(qryBlock, workerNum);
			V **qryKHistList = get<0>(qryTuble);
			uint64_t **qryMonoHistList = get<1>(qryTuble);
			std::string **qryInfoList = get<2>(qryTuble);
			int *qryLenList = get<3>(qryTuble);

			for (int i = 0; i < dbSize; i++) {
				std::string dbSeq = *dbInfoList[i];
				double *v = id->score(dbKHistList[i], qryKHistList,
						dbMonoHistList[i], qryMonoHistList, qrySize, workerNum,
						dbLenList[i], qryLenList);

				for (int j = 0; j < qrySize; j++) {
					if (v[j] > 0.0) {
						out << dbSeq << "\t" << *qryInfoList[j] << "\t" << v[j]
								<< std::endl;
					}
				}
				delete[] v;
			}
			// Free query block
			id->freeBlock(qryTuble, qrySize, workerNum);
		}
		// Free database block
		id->freeBlock(dbTuble, dbSize, workerNum);
	}
	cout << endl;

	// Close output file.
	out.flush();
	out.close();
	delete id;
}

template<class V>
void ReaderAlignerCoordinator::helper2(string fileDb, string fileQry,
		string fileOut, string dlm, bool isAllVsAll, DataGenerator *g,
		Serializer *serializer) {

	bool canSkip = !canReportAll;

	IdentityCalculator<V> *id;
	if (canFillModel) {
		id = new IdentityCalculator<V>(*serializer, threshold, canSkip,
				canRelax);
	} else if (canSaveModel) {
		id = new IdentityCalculator<V>(g, workerNum, threshold, canSkip,
				canRelax, modelFile);
	} else {
		id = new IdentityCalculator<V>(g, workerNum, threshold, canSkip,
				canRelax);
	}

	/**
	 * If the output file is not provided, identity is used in the training mode only
	 */
	if (fileOut.empty()) {
		std::cout << "Training is done." << std::endl;
		return;
	}

	if (canRelax) {
		std::cout << "Relaxing the threshold" << std::endl;
	}

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
		vector<Aligner<V>*> alignerList;
		alignerList.reserve(workerNum);
		vector<future<std::pair<bool, stringstream*> > > futureList;
		futureList.reserve(workerNum);
		for (int i = 0; i < workerNum; i++) {
			Aligner<V> *aligner = new Aligner<V>(*id, dbBlock, dlm,
					canReportAll, threshold, canRelax);
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
	out.flush();
	out.close();
	delete id;
}
