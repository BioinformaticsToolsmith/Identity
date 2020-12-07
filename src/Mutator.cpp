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
 * Mutator.cpp
 *
 *  Created on: Oct 23, 2019
 *      Author: Hani Zakaria Girgis, PhD
 *
 */

#include "Mutator.h"
#include <iostream>

/**
 * The composition list is constructed from the input sequence
 */
Mutator::Mutator(const string *oSequenceIn, int maxBlockIn, int seed,
		int minBlockIn) {
	oSequence = oSequenceIn;
	makeCompositionList();
	help(maxBlockIn, seed, minBlockIn);
}

/**
 * The composition list is provided by the client
 */
Mutator::Mutator(const string *oSequenceIn, int maxBlockIn, int seed,
		vector<double> *compositionListIn, int minBlockIn) {
	oSequence = oSequenceIn;
	compositionList = compositionListIn;
	ownCompositionList = false;
	help(maxBlockIn, seed, minBlockIn);
}

void Mutator::help(int maxBlockIn, int seed, int minBlockIn) {
	// Preconditions
	if (maxBlockIn <= 1) {
		cerr << "The maximum block size must be greater than 1" << endl;
		cerr << "Maximum size: " << maxBlockIn << endl;
		throw std::exception();
	}

	if (minBlockIn <= 1) {
		cerr << "The minimum block size must be greater than 1" << endl;
		cerr << "Minimum size: " << minBlockIn << endl;
		throw std::exception();
	}

	if (minBlockIn > maxBlockIn) {
		cerr << "The minimum block size cannot be greater than its maximum."
				<< endl;
		cerr << "Minimum size: " << minBlockIn << ", maximum size: "
				<< maxBlockIn << endl;
		throw std::exception();
	}
	// End pre-conditions

	inversionFactor = Parameters::getInversionFactor();
	translocationFactor = Parameters::getTranslocationFactor();

	// Segment coordinates are inclusive [s,e]
	// This peace of code is duplicated in KmerHistogram
	unknown = Parameters::getUnknown();
	segmentList = new vector<pair<int, int> >();
	// Store non-N/non-X index
	int start = -1;
	int oLen = oSequence->size();
	for (unsigned int i = 0; i < oLen; i++) {
		if (oSequence->at(i) != unknown && start == -1) {
			start = i;
		} else if (oSequence->at(i) == unknown && start != -1) {
			segmentList->push_back(make_pair(start, i - 1));
			start = -1;
		} else if (i == oLen - 1 && oSequence->at(i) != unknown
				&& start != -1) {
			segmentList->push_back(make_pair(start, i));
			start = -1;
		}
	}
	// Post condition
	if (segmentList->empty()) {
		cerr << "At least one valid segment is required." << endl;
		throw std::exception();
	}
	// Effective length is the number of valid nucleotides or a.a.
	for (auto p : *segmentList) {
		// cout << "Segment: " << p.first << " " << p.second << endl;
		effectiveLength += p.second - p.first + 1;
	}

	minBlock = minBlockIn;
	maxBlock = maxBlockIn;
	if (maxBlock < minBlock) {
		maxBlock = minBlock;
	}

	mutationList = new vector<int>();
	mutationList->reserve(Mutation::TRANSLOCATION);

	// The following lines are needed for generating random nucleotide
	aLimit = compositionList->at(0);
	cLimit = compositionList->at(0) + compositionList->at(1);
	gLimit = compositionList->at(0) + compositionList->at(1)
			+ compositionList->at(2);
	tLimit = 1.0;

	// Credit: https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	gen = new std::default_random_engine(seed);
	zeroOneRand = new uniform_real_distribution<>(0.0, 1.0);
}

void Mutator::makeCompositionList() {
	// Calculate the nucleotide composition of the original sequence
	compositionList = new vector<double>(4, 0.0);
	ownCompositionList = true;
	int oSeqSize = oSequence->size();
	// No need to check the lower case because the reader converts all
	// nucleotide to large case.
	for (int i = 0; i < oSeqSize; i++) {
		switch (oSequence->at(i)) {
		case 'A':
		case 'a':
			compositionList->at(0)++;break
;			case 'C':
			case 'c':
			compositionList->at(1)++; break;
			case 'G':
			case 'g':
			compositionList->at(2)++; break;
			case 'T':
			case 't':
			compositionList->at(3)++; break;
		}
	}
	for (int i = 0; i < 4; i++) {
		compositionList->at(i) /= oSeqSize;
	}
}

void Mutator::enableSinglePoint() {
	//mutationList->push_back(Mutation::INSERTION);
	mutationList->push_back(Mutation::DELETION);
	mutationList->push_back(Mutation::MISMATCH);
}

void Mutator::enableBlock() {
	//mutationList->push_back(Mutation::B_INSERTION);
	mutationList->push_back(Mutation::B_DELETION);
	mutationList->push_back(Mutation::DUPLICATION);
}

void Mutator::enableInverstion() {
	mutationList->push_back(Mutation::INVERSION);
}

void Mutator::enableTranslocation() {
	mutationList->push_back(Mutation::TRANSLOCATION);
	// If no block deletion is already in the list, add it
	if (std::count(mutationList->begin(), mutationList->end(),
			Mutation::B_DELETION) == 0) {
		mutationList->push_back(Mutation::B_DELETION);
	}
}

Mutator::~Mutator() {
	mutationList->clear();
	delete mutationList;

	if (ownCompositionList) {
		compositionList->clear();
		delete compositionList;
	}

	segmentList->clear();
	delete segmentList;

	delete gen;
	delete zeroOneRand;
}

/*
 * Return a nucleotide according to the distribution of the original sequence.
 */
char Mutator::getRandomNucleotide() {
	double p = (*zeroOneRand)(*gen);
	char r;
	if (p <= aLimit) {
		r = 'A';
	} else if (p <= cLimit) {
		r = 'C';
	} else if (p <= gLimit) {
		r = 'G';
	} else {
		r = 'T';
	}
	return r;
}

/**
 * Thread unsafe
 *
 * Algorithm outline:
 * . The algorithm is linear.
 * . Single point mutation types: insertion, deletion, and mismatch.
 * . Block mutation types: insertion, deletion, duplication, translocation, and inversion.
 * . Allowing a mix of these two types while restricting the size of the block to few nucleotides can
 * 		produce a biologically relevant mutation model.
 * . Block size is selected randomly from a range provided by the client.
 * . The client must provide the maximum block size.
 * . The minimum block size has a default value of 2.
 * . These mutation are not overlapping.
 * . Mutations are generated based on the nucleotide composition of the input sequence or provided by the client.
 * . How to construct the mutated string:
 *		. Determine the number of nucleotides to be mutated.
 *		. Distribute these nucleotides uniformly over the input sequence. Call these positions indexes.
 *		. If a block mutation of size x, then move x indexes, i.e. the next mutation happens after x indexes.
 * . Translocation is treaded as deletion followed by insertion. The algorithm saves the coordinates of each deleted segment.
 * 		Then when a translocation is chosen, the algorithm selects the last deleted segment and inserts it at the new index.
 * 		If no deleted segments are available, the algorithm chooses another mutation type.
 * . Randomization has been fixed, i.e. the same sequences are produced using the same seed.
 */
pair<string*, double> Mutator::mutateSequence(double mutationRate) {
	// Pre-conditions
	if (mutationRate < 0 || mutationRate > 1) {
		cerr << "Mutation rate must be >= 0 and <= 1, but received ";
		cerr << mutationRate << endl;
		throw std::exception();
	}

	if (mutationList->empty()) {
		cerr << "No mutation types were specified" << endl;
		throw std::exception();
	}
	// End of pre-conditions

	int mutationTotal = round(mutationRate * oSequence->length());
	if (mutationTotal > effectiveLength) {
		mutationTotal = effectiveLength;
	}

	int oLength = oSequence->length();

	string *mSequence = new string("");

	// Delete if good
	// mSequenceList->push_back(mSequence);

	// Extra space is trimmed at the end of this method.
	mSequence->reserve(oSequence->size() + mutationTotal);

	// This vector holds the locations of the deleted segments.
	vector<pair<int, int> > deleteList;
	// Initially indexes to be mutated are distributed uniformly.
	// If single point mutations are enabled only, the resulting
	// mutations are uniformly distributed. As a remedy, the client
	// should enable block mutations with a small maximum length, e.g. 5.
	double interval = effectiveLength / (double) mutationTotal;
	// cout << "Interval: " << interval << endl;

	int mutationNum = mutationList->size();
	// Initialize random number generators
	// Random number generator for mutation type
	std::uniform_int_distribution<int> randMutIndexDist(0, mutationNum - 1);
	// Random number generator for block size
	std::uniform_int_distribution<int> randBlockSizeDist(minBlock, maxBlock);

	int mutationRemaining = mutationTotal;

	double alignLen = oLength;
	double matchNum = oLength;

	double identity = 0.0;
	if (mutationTotal < 1) {
		mSequence->append(*oSequence);
		identity = 1.0;
	} else {
		int segIndex = 0;
		pair<int, int> seg = segmentList->at(segIndex);
		int segStart = seg.first;
		int segEnd = seg.second;
		int segNum = segmentList->size();

		int skipped = 0;
		for (int i = 0; i - skipped < mutationTotal;) {
			int index = i * interval;
			// It got to the end of a sequence, but it did not use up all available mutations.
			if (index >= oLength) {
				break;
			}

			int oldIndex = index;

			// Check and skip indexes that are in an invalid region
			while (i < mutationTotal && index < oLength
					&& oSequence->at(index) == unknown) {
				index = ++i * interval;
				skipped++;
			}
			// If no valid index found
			if (index >= oLength || oSequence->at(index) == unknown) {
				mSequence->append(
						oSequence->substr(oldIndex, oLength - oldIndex));
				break;
			}
			// If some indexes are skipped, copy the in-between region
			if (index != oldIndex) {
				mSequence->append(
						oSequence->substr(oldIndex, index - oldIndex));
			}
			// Find the correct segment. The end is needed to limit the
			// random block size.
			while (!(index >= segStart && index <= segEnd)
					&& segIndex < segNum - 1) {
				pair<int, int> seg = segmentList->at(++segIndex);
				segStart = seg.first;
				segEnd = seg.second;
			}
			// Test condition
			if (!(index >= segStart && index <= segEnd)) {
				cerr << "Index not in segment" << endl;
				throw std::exception();
			}

			// Select mutation type randomly
			int mutationType = mutationList->at(randMutIndexDist(*gen));
			// Translocation depends on block deletion.
			// Choose another mutation type if translocation is selected but
			// the delete list is empty.
			while (mutationType == Mutation::TRANSLOCATION && deleteList.empty()) {
				mutationType = mutationList->at(randMutIndexDist(*gen));
			}

			int randBlockSize;
			int nextIndex;
			// Variables for block mutations
			if (mutationType > Mutation::MISMATCH) {
				// Calculate random block size minBlock–-maxBlock
				randBlockSize = randBlockSizeDist(*gen);

				// Make sure that the random block size is not greater than the number
				// of nucleotides need to be mutated.
				if (randBlockSize > mutationRemaining) {
					randBlockSize = mutationRemaining;
				}
				// Make sure that the random block size does not go beyond the end of the segment.
				if (index + randBlockSize > segEnd + 1) {
					randBlockSize = segEnd - index + 1;
				}

				nextIndex = (i + randBlockSize) * interval;
				mutationRemaining -= randBlockSize;
			}
			// Variables for single point mutations
			else {
				randBlockSize = 1;
				nextIndex = (i + 1) * interval;
				mutationRemaining--;
			}
			// If this is the last mutation, then the next index is the size
			// of the original string exclusive. This variable (nextIndex) is
			// exclusive.
			if (mutationRemaining == 0 || nextIndex > oLength) {
				nextIndex = oLength;
			}

			if (nextIndex < index) {
				std::cerr
						<< "Error in Mutator: The next index cannot be less than the current index.";
				std::cerr << std::endl;
				throw std::exception();
			}

			// Skip as many indexes as the block size
			bool isBlockMutation = true;
			switch (mutationType) {
			case Mutation::INSERTION: {
				isBlockMutation = false;
				// Insert a random nucleotide
				mSequence->append(1, getRandomNucleotide());
				// Copy from index to the next index - 1
				mSequence->append(oSequence->substr(index, nextIndex - index));
				// Update alignment length only
				alignLen += 1;
			}
				break;
			case Mutation::DELETION: {
				isBlockMutation = false;

				// Copy from index+1 to the next index - 1
				mSequence->append(
						oSequence->substr(index + 1, nextIndex - index - 1));
				// Update number of matches only
				matchNum -= 1;
			}
				break;
			case Mutation::MISMATCH: {
				isBlockMutation = false;
				char randChar = getRandomNucleotide();
				mSequence->append(1, randChar);
				// Copy from index to the next index - 1
				mSequence->append(
						oSequence->substr(index + 1, nextIndex - index - 1));

				// ToDo: check to see if the new character is the same as
				// the old one and update the number of matches. I allowed this
				// to skip a mutation spot.

				// Update the number of matches only
				if (randChar != oSequence->at(index)) {
					matchNum -= 1;
				}
			}
				break;
			case Mutation::B_INSERTION: {
				for (int h = 0; h < randBlockSize; h++) {
					mSequence->append(1, getRandomNucleotide());
				}

				mSequence->append(oSequence->substr(index, nextIndex - index));
				// Update alignment length only
				alignLen += randBlockSize;
			}
				break;
			case Mutation::B_DELETION: {
				mSequence->append(
						oSequence->substr(index + randBlockSize,
								nextIndex - index - randBlockSize));
				deleteList.push_back(
						make_pair(index, index + randBlockSize - 1));
				// Update number of matches only
				matchNum -= randBlockSize;
			}
				break;
			case Mutation::DUPLICATION: {
				// Copy the block
				mSequence->append(oSequence->substr(index, randBlockSize));

				// Copy the segment starting with the block
				mSequence->append(oSequence->substr(index, nextIndex - index));
				// Update alignment length only
				alignLen += randBlockSize;
			}
				break;
			case Mutation::INVERSION: {
				// Invert the segment
				string seg("");
				seg.reserve(randBlockSize);
				for (int h = index + randBlockSize - 1; h >= index; h--) {
					seg.append(1, oSequence->at(h));
				}

				mSequence->append(seg);

				mSequence->append(
						oSequence->substr(index + randBlockSize,
								nextIndex - index - randBlockSize));

				// Inversion should be treated in a similar way to mismatch.
				// However, it should consider that fact that all nucleotides
				// are still present, but in the opposite direction.
				double misMatch = 0.0; // Number of mismatches
				for (int h = index; h < index + randBlockSize; h++) {
					if (oSequence->at(h) != seg.at(h - index)) {
						misMatch++;
					}
				}
				// Decrement by 50% of the mismatches
				matchNum -= round(misMatch * inversionFactor);
			}
				break;
			case Mutation::TRANSLOCATION: {
				// Note: a tranlocation undoes a previous deletion
				pair<int, int> segment = deleteList.back();
				// Copy a deleted block

				mSequence->append(
						oSequence->substr(segment.first,
								segment.second - segment.first + 1));

				// Copy the original segment
				mSequence->append(oSequence->substr(index, nextIndex - index));
				deleteList.pop_back();

				// Update alignment length and match number
				// The final effect of a translocation (deletion followed by insertion)
				// should have no effect on the alignment length. Recall, deleting a block
				// does not affect this variable. So, nothing to be done. However, deleting
				// a block reduces the number of matches. Therefore, this effect needs to
				// be undone. One way is to undo it completely, but that would be unrealistic.
				// Another arbitrary yet reasonable choice is to undo it by 50%.
				matchNum += round(randBlockSize * translocationFactor);
			}
				break;
			default: {
				cerr << "Undefined mutation type" << endl;
				throw std::exception();
			}
			}

			// Update the next index and the maximum allowed block size.
			if (isBlockMutation) {
				i += randBlockSize;
			} else {
				i++;
			}
		}
		identity = matchNum / alignLen;
	}
	// Trim extra space that was allocated at the beginning of this method
	mSequence->shrink_to_fit();
	// Convert the mutated string to upper case
	// Based on code found at https://thispointer.com
	std::for_each(mSequence->begin(), mSequence->end(), [](char &c) {
		if (c >= 97)
			c = ::toupper(c);
	});

	return std::make_pair(mSequence, identity);
}
