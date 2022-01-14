/*
 Identity calculates DNA sequence identity scores rapidly without alignment.

 Copyright (C) 2020-2022 Hani Z. Girgis, PhD

 Academic use: Affero General Public License version 1.

 Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 */

/*
 * Reservoir.cpp
 *
 *  Created on: Jul 5, 2021
 *      Author: Hani Z. Girgis, Phd
 */

template<class V>
Reservoir<V>::Reservoir() {
	seed = 17;
}

template<class V>
Reservoir<V>::~Reservoir() {
	if (!kHistList.empty() || !monoHistList.empty() || !infoList.empty()
			|| !lenList.empty()) {
		cerr << "Warning: The reservoir is destroyed, but it is not empty!";
		cerr << endl;
	}
}

template<class V>
void Reservoir<V>::add(
		std::tuple<V**, uint64_t**, std::string**, int*, int> t) {
	V **l1 = get < 0 > (t);
	uint64_t **l2 = get < 1 > (t);
	std::string **l3 = get < 2 > (t);
	int *l4 = get < 3 > (t);
	int blockSize = get < 4 > (t);

	for (int i = 0; i < blockSize; i++) {
		kHistList.push_back(l1[i]);
		monoHistList.push_back(l2[i]);
		infoList.push_back(l3[i]);
		lenList.push_back(l4[i]);
	}

	delete[] l1;
	delete[] l2;
	delete[] l3;
	delete[] l4;
}

template<class V>
std::tuple<V**, uint64_t**, std::string**, int*, int> Reservoir<V>::remove(
		int blockSize) {

	shuffle();

	if (blockSize > kHistList.size()) {
		blockSize = kHistList.size();
	}
	V **l1 = new V*[blockSize];
	uint64_t **l2 = new uint64_t*[blockSize];
	std::string **l3 = new std::string*[blockSize];
	int *l4 = new int[blockSize];

	for (int i = 0; i < blockSize; i++) {
		l1[i] = kHistList[i];
		l2[i] = monoHistList[i];
		l3[i] = infoList[i];
		l4[i] = lenList[i];
	}

	kHistList.erase(kHistList.begin(), kHistList.begin() + blockSize);
	monoHistList.erase(monoHistList.begin(), monoHistList.begin() + blockSize);
	infoList.erase(infoList.begin(), infoList.begin() + blockSize);
	lenList.erase(lenList.begin(), lenList.begin() + blockSize);

	return std::make_tuple(l1, l2, l3, l4, blockSize);
}

template<class V>
void Reservoir<V>::shuffle() {
	int s = size();

	// Make a list of indexes
	vector<int> indexList(s, -1);
	for (int i = 0; i < s; i++) {
		indexList[i] = i;
	}

	// Shuffle the index list
	std::shuffle(indexList.begin(), indexList.end(),
			std::default_random_engine(seed));
	seed++;

	// Make new vectors
	std::vector<V*> kHistListShuffled(s, 0);
	std::vector<uint64_t*> monoHistListShuffled(s, 0);
	std::vector<std::string*> infoListShuffled(s, nullptr);
	std::vector<int> lenListShuffled(s, 0);

	// Copy items according to the shuffled indexes
	for (int i = 0; i < s; i++) {
		kHistListShuffled[i] = kHistList[indexList[i]];
		monoHistListShuffled[i] = monoHistList[indexList[i]];
		infoListShuffled[i] = infoList[indexList[i]];
		lenListShuffled[i] = lenList[indexList[i]];
	}

	// Copy the shuffled vectors to the old ones
	kHistList.clear();
	monoHistList.clear();
	infoList.clear();
	lenList.clear();
	kHistList = kHistListShuffled;
	monoHistList = monoHistListShuffled;
	infoList = infoListShuffled;
	lenList = lenListShuffled;
}

template<class V>
int Reservoir<V>::size() {
	return kHistList.size();
}
