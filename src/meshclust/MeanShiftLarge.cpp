/*
 MeShClust v3.0 clusters sequences using the mean shift algorithm and alignment-free identity scores.

 Copyright (C) 2020-2022 Hani Z. Girgis, PhD

 Academic use: Affero General Public License version 1.

 Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 */

/*
 * MeanShiftLarge.cpp
 *
 *  Created on: Dec 23, 2020
 *      Author: Hani Zakaria Girgis, PhD
 */
template<class V>
MeanShiftLarge<V>::MeanShiftLarge(string file, int size, int v, int p,
		IdentityCalculator<V> &id, int t, bool yesNo, string out, double th,
		bool canEval, bool relax) :
		dbFile(file), blockSize(size), vBlockSize(v), passNum(p), identity(id), threadNum(
				t), canAssignAll(yesNo), outFile(out), threshold(th), canEvaluate(
				canEval), canRelax(relax) {

	seqNum = 0;

	clusterReservoir();
}

template<class V>
MeanShiftLarge<V>::~MeanShiftLarge() {
	delete ms;
}

template<class V>
void MeanShiftLarge<V>::clusterReservoir() {
	std::cout << std::endl << "Clustering ... " << std::endl;
	FastaReader reader(dbFile, blockSize);
	auto block = reader.read();
	seqNum += block->size();

	std::cout << std::endl << "Data run 1 ..." << std::endl;
	ms = new MeanShift<V>(block, identity, threadNum, threshold);
	ms->removeSingles();
	// printStatus(true);
	printStatus();

	if (reader.isStillReading()) {
		ms->updateAccumulatedMean();
		MeanShift < V > *cSinglesMs = nullptr;
		MeanShift < V > *pSinglesMs = nullptr;
		bool canAddCenters = false;
		reservoir.add(ms->findUnassignedData());

		for (int i = 0; i < passNum; i++) {
			int clustNum = ms->getClusterList()->size();

			if (i > 0) {
				reader.setBlockSize(vBlockSize);
				// std::cout << "New block size is: " << vBlockSize << std::endl;
				std::cout << std::endl << "Data run " << (i + 1) << " ..."
						<< std::endl;
			}

			bool isReading = reader.isStillReading();
			bool isFull = reservoir.size() > 0; // Initially false? No

			while (isReading || isFull) {
				if (isReading) {
					ms->clearData();

					auto block = reader.read();
					int size = block->size();
					seqNum += size;

					ms->updateReferenceData(block);
					isReading = reader.isStillReading();
				}

				if (canAddCenters) {
					ms->addClusters(pSinglesMs->getClusterList());
				}

				if (isReading || canAddCenters) {
					ms->run(1, i == 0 ? true : false);
					ms->updateAccumulatedMean();
					if (i == 0) {
						auto u = ms->findUnassignedData();
						reservoir.add(u);
						// Determine the next block size for reading sequences
						int n = vBlockSize;
						if (get < 4 > (u) > 0) {
							n = blockSize * blockSize / (double) get < 4 > (u);
						}
						reader.setBlockSize(n > vBlockSize ? vBlockSize : n);
					}
					// printStatus(true);
					printStatus();
				}

				if (reservoir.size() > blockSize
						|| (!isReading && reservoir.size() > 0)) {
					cSinglesMs = new MeanShift<V>(reservoir.remove(blockSize),
							identity, threadNum, threshold);
					cSinglesMs->removeSingles();
					if (isReading) {
						reservoir.add(cSinglesMs->findUnassignedData());
					}
					if (pSinglesMs != nullptr) {
						delete pSinglesMs;
					}
					pSinglesMs = cSinglesMs;
					canAddCenters = true;
				} else {
					canAddCenters = false;
				}

				// Update conditions
				isReading = reader.isStillReading();
				isFull = reservoir.size() > 0;
				//printStatus();
			}

			if (canAddCenters) {
				ms->addClusters(pSinglesMs->getClusterList());
				ms->run(1, false);
				ms->updateAccumulatedMean();

				// printStatus(true);
				printStatus();
			}

			// Post-condition
			if (reservoir.size() > 0) {
				std::cerr << "The reservoir must be empty, but its size is: ";
				std::cerr << reservoir.size() << std::endl;
				throw std::exception();
			}

			printStatus(true);

			// Free memory
			if (pSinglesMs != nullptr) {
				delete pSinglesMs;
				pSinglesMs = nullptr;
			}

			// Restart
			reader.restart();
			seqNum = 0;
			printCounter = 0;
			canAddCenters = false;

			if (i > 0 && clustNum == ms->getClusterList()->size()) {
				break;
			}
		}
	}
	ms->clearData();

	// Assign data
	assign();
}

template<class V>
void MeanShiftLarge<V>::printStatus(bool printAnyWay) {
	if (printAnyWay
			|| seqNum >= (printCounter * Parameters::getMsPrintBlock())) { // 50000
		if (!printAnyWay) {
			printCounter++;
		}
		std::cout << "\tProcessed sequences: " << seqNum << std::endl;
		std::cout << "\tUnprocessed sequences: " << reservoir.size()
				<< std::endl;
		std::cout << "\tFound centers: " << ms->getClusterList()->size()
				<< std::endl;
	}
}

template<class V>
Matrix MeanShiftLarge<V>::calcAllCenterVsAllCenter() {
	return ms->calcAllCenterVsAllCenter();
}

/**
 * Assign each sequence in the database to a cluster
 * Write results
 */
template<class V>
void MeanShiftLarge<V>::assign() {
	std::cout << std::endl << "Assigning ..." << std::endl;

	double error = identity.getError();
	double relaxThreshold = threshold - error;

	// All identities are calculated regardless of the threshold.
	if (canEvaluate || canAssignAll) {
		identity.setCanSkip(false);
	}

	// Work on the original file
	FastaReader reader(dbFile, vBlockSize);

	auto clusterList = ms->getClusterList();
	int clusterNum = clusterList->size();

	// This call instantiate the member list
	vector<ClusterInfo*> clusterInfoList;
	int identifier = 0;
	clusterInfoList.reserve(clusterNum);

	for (int i = 0; i < clusterNum; i++) {
		clusterInfoList.push_back(new ClusterInfo(++identifier));
	}

	// Holds singles
	vector<ClusterInfo*> singleList;

	int dataSize = 0;
	while (reader.isStillReading()) {
		Block *block = reader.read();
		int blockSize = block->size();
		dataSize += blockSize;
		double *res[clusterNum];
		auto tup = identity.unpackBlock(block, threadNum);
		auto kHistList = std::get < 0 > (tup);
		auto monoHistList = std::get < 1 > (tup);
		auto infoList = std::get < 2 > (tup);
		auto lenList = std::get < 3 > (tup);

		for (int i = 0; i < clusterNum; i++) {
			auto c = clusterList->at(i);
			res[i] = identity.score(c->getKHistMean(), kHistList,
					c->getMonoHistMean(), monoHistList, blockSize, threadNum,
					c->getLength(), lenList);
		}

		typedef struct {
			int index;
			double scoreWithCenter;
			double scoreWithNeighbor;
			int membership;
		} info;

		info assignment[blockSize];
#pragma omp parallel for schedule(static) num_threads(threadNum)
		for (int i = 0; i < blockSize; i++) {
			double max = -1.0;
			int index = -1;
			for (int j = 0; j < clusterNum; j++) {
				double val = res[j][i];
				if (val > max) {
					max = val;
					index = j;
				}
			}

			// Fill scores with the closest cluster
			double secondBest = -1.0;
			if (canEvaluate) {
				for (int j = 0; j < clusterNum; j++) {
					double val = res[j][i];
					if (val < max && val > secondBest) {
						secondBest = val;
					}
				}
			}

			if (max >= threshold) {
				assignment[i] =
						{ index, max, secondBest, ClusteringUtil::MEMBER };
			} else {
				if (canRelax && max >= relaxThreshold) {
					assignment[i] = { index, max, secondBest,
							ClusteringUtil::EXTENDED };
				} else if (canAssignAll) {
					assignment[i] = { index, max, secondBest,
							ClusteringUtil::OUTSIDE };
				} else {
					// This is a singleton
					assignment[i] = { -1, 1.0000, max, ClusteringUtil::MEMBER };
				}
			}
		}

		for (int i = 0; i < blockSize; i++) {
			int index = assignment[i].index;
			if (index >= 0) {
				clusterInfoList[index]->addMember(infoList[i],
						assignment[i].scoreWithCenter,
						assignment[i].scoreWithNeighbor,
						assignment[i].membership);
			} else {
				auto clusterInfo = new ClusterInfo(++identifier);
				clusterInfo->addMember(infoList[i],
						assignment[i].scoreWithCenter,
						assignment[i].scoreWithNeighbor,
						assignment[i].membership);
				singleList.push_back(clusterInfo);
			}
		}

		// Free memory
#pragma omp parallel for schedule(static) num_threads(threadNum)
		for (int i = 0; i < clusterNum; i++) {
			delete[] res[i];
		}

#pragma omp parallel for schedule(static) num_threads(threadNum)
		for (int i = 0; i < blockSize; i++) {
			delete[] kHistList[i];
			delete[] monoHistList[i];
			// Do not delete the info (header) because it is still is use.
			// It will be deleted in ClusterInfo
		}

		delete[] kHistList;
		delete[] monoHistList;
		delete[] infoList;
		delete[] lenList;

		if ((assignCounter * Parameters::getMsPrintBlock()) <= dataSize) {
			std::cout << "\tSequences assigned to clusters: " << dataSize;
			std::cout << std::endl;
			assignCounter++;
		}
	}

	// Remove empty clusters
	clusterInfoList.erase(
			std::remove_if(clusterInfoList.begin(), clusterInfoList.end(),
					[](ClusterInfo *ptr) {
						return ptr->getSize() == 0;
					}), clusterInfoList.end());
	int newSize = clusterInfoList.size();

	// Reset the identifiers
	if (clusterNum > newSize) {
		for (int u = 0; u < newSize; u++) {
			clusterInfoList[u]->setIdentifier(u + 1);
		}
	}

	if (canEvaluate) {
		std::cout << std::endl << "Evaluating ..." << std::endl;
		Matrix ava = calcAllCenterVsAllCenter();
		ClusterEvaluator evaluator(ava, clusterInfoList, dataSize);
		evaluator.printAll(true);
	}

	// Write out results
	ofstream out(outFile);
	for (int i = 0; i < newSize; i++) {
		out << clusterInfoList[i]->toString() << endl;
	}

	int singleCount = singleList.size();
	for (int j = 0; j < singleCount; j++) {
		out << singleList.at(j)->toString() << endl;
		delete singleList.at(j);
	}
	out.close();

#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < newSize; i++) {
		delete clusterInfoList.at(i);
	}
}
