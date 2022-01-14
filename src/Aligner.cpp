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
 * Aligner.cpp
 *
 *  Created on: Nov 15, 2019
 *      Author: Hani Zakaria Girgis, PhD
 *  An instance of this class calculates the All-vs-All on a block of sequences.
 *  It is thread safe designed for parallel execution.
 *
 */

/**
 * Block a will NOT be deleted here because it is
 * being processed by other threads as well.
 */
template<class V>
Aligner<V>::Aligner(IdentityCalculator<V> &c, Block *a, string dlmIn,
		bool filter, double cutoff, bool canRelax) :
		identity(c) {
	blockA = a;
	dlm = dlmIn;
	threshold = cutoff;

	canReportAll = filter;
	ssPtr = new stringstream();

	if (canRelax) {
		error = identity.getError();
	}

	k = identity.getK();
}

template<class V>
Aligner<V>::~Aligner() {
	delete ssPtr;

	if (buffer.size() > 0) {
		std::cerr << "Aligner error: Queue must be empty. " << std::endl;
		std::cerr << "Queue size is: " << buffer.size() << std::endl;
	}
}

template<class V>
pair<bool, stringstream*> Aligner<V>::start() {
	// Keep processing blocks as they are enqueued.
	while (true) {
		if (buffer.size() > 0) {
			processBlock();
		} else if (canStop && buffer.size() == 0) {
			break;
		}
	}

	return std::make_pair(canWrite, ssPtr);
}

template<class V>
pair<bool, stringstream*> Aligner<V>::getResults() {
	// Keep processing blocks as they are enqueued.
	return std::make_pair(canWrite, ssPtr);
}

/**
 * Thread safe
 * Note: This block and its contents will be deleted
 * after processing it.
 */
template<class V>
void Aligner<V>::enqueueBlock(pair<Block*, bool> p) {
	buffer.push(p);
}

/**
 * Thread safe
 * Block A: Query
 * Block B: Database
 */
template<class V>
void Aligner<V>::processBlock() {
	// Get the front of the the queue, but do not pop it yet.
	auto blockB = buffer.front().first;
	int sizeA = blockA->size();
	int sizeB = blockB->size();

	static KmerHistogram<uint64_t, V> kTable(k);
	static KmerHistogram<uint64_t, uint64_t> monoTable(1);

	for (int j = 0; j < sizeA; j++) {
		int init = 0;
		// If the two blocks have the same contents.
		if (buffer.front().second) {
			init = j + 1;
		}
		auto p1 = blockA->at(j);

		string *info1 = p1.first;
		string *seq1 = p1.second;

		V *h1 = kTable.build(seq1);
		uint64_t *mono1 = monoTable.build(seq1);

		double l1 = seq1->size();
		int sizeB = blockB->size();

		for (int hani = init; hani < sizeB; hani++) {
			auto p2 = blockB->at(hani);
			string *seq2 = p2.second;
			int l2 = seq2->size();

			double ratio = l1 < l2 ? l1 / l2 : l2 / l1;
			if (!canReportAll && ratio < threshold) {
				continue;
			}

			V *h2 = kTable.build(seq2);
			uint64_t *mono2 = monoTable.build(seq2);

			double res = identity.score(h1, h2, mono1, mono2, ratio, l1, l2);

			if (canReportAll || res > 0.0) {
				canWrite = true;
				(*ssPtr) << *info1 << dlm << *p2.first << dlm
						<< std::setprecision(4) << res << std::endl;
			}

			delete[] h2;
			delete[] mono2;
		}

		delete[] h1;
		delete[] mono1;
	}

	// Pop the block and free its memory
	FastaReader::deleteBlock(blockB);
	buffer.pop();
}

template<class V>
int Aligner<V>::getQueueSize() {
	return buffer.size();
}

/**
 * Thread safe
 */
template<class V>
void Aligner<V>::stop() {
	canStop = true;
}
