/*
 * KmerHistogram.cpp
 *
 * Created on: Jul 25, 2012
 * Modified on: 12/21/2019
 * Author: Hani Zakaria Girgis, PhD
 */
template<class I, class V>
KmerHistogram<I, V>::KmerHistogram(int keyLength) :
		k(keyLength) {
	if (k < 1) {
		cerr << "k must be at least 1." << endl;
		throw std::exception();
	}

	// Initialize digit list
	for (int i = 0; i <= 'T'; i++) {
		digitList[i] = '*';
	}
	digitList['C'] = 0;
	digitList['T'] = 1;
	digitList['A'] = 2;
	digitList['G'] = 3;

	// Initialize bases
	// bases = new I[k];
	for (int i = k - 1; i >= 0; i--) {
		bases[k - 1 - i] = (I) pow(4.0, i);
	}

	// Initialize mMinusOne
	// mMinusOne = new I[4];
	for (int i = 0; i < 4; i++) {
		mMinusOne[i] = i * bases[0];
	}

	// Get maximum size of table
	string temp;
	for (int i = 0; i < k; i++) {
		temp.append(1, 'G');
	}
	maxTableSize = hash(&temp) + 1;
}

template<class I, class V>
KmerHistogram<I, V>::~KmerHistogram() {
	// delete[] bases;
	// delete[] mMinusOne;
}

/**
 * key: The first k nucleotides comprise the key.
 */
template<class I, class V>
I KmerHistogram<I, V>::hash(const string *key) {
	return hash(key, 0);
}

/**
 * sequence: An array of valid nucleotides.
 * keyStart: The start index of the key.
 */
template<class I, class V>
I KmerHistogram<I, V>::hash(const string *sequence, int keyStart) {
	I index = 0;
	const char *arr = sequence->c_str();
	for (int i = 0; i < k; i++) {
		index += bases[i] * digitList[arr[keyStart + i]];
	}

	return index;
}

/**
 *
 */
template<class I, class V>
void KmerHistogram<I, V>::hash(const string *sequence, int start, int end,
		vector<I> *hashList) {

	I lastHash = hash(sequence, start);
	hashList->push_back(lastHash);
	const char *arr = sequence->c_str();
	for (int i = start + 1; i <= end; i++) {
		I s1 = 4 * (lastHash - mMinusOne[digitList[arr[i - 1]]])
				+ digitList[arr[i + k - 1]];
		hashList->push_back(s1);
		lastHash = s1;
	}
}

/**
 * Build a k-mer histogram from DNA sequence efficiently.
 * sequence: Count k-mers in this sequence.
 * Memory: This method builds a histogram and puts it on the heap.
 * 	The client is responsible for destroying the histogram.
 */
template<class I, class V>
V* KmerHistogram<I, V>::build(const string *sequence) {
	// Segment coordinates are inclusive [s,e]
	// This peace of code is duplicated in KmerHistogram
	char unknown = Parameters::getUnknown();
	vector < pair<int, int> > segmentList;
	// Store non-N/non-X index
	int start = -1;
	int len = sequence->size();
	const char *array = sequence->c_str();
	for (unsigned int i = 0; i < len; i++) {
		char c = array[i];
		if (c != unknown && start == -1) {
			start = i;
		} else if (c == unknown && start != -1) {
			segmentList.push_back(make_pair(start, i - 1));
			start = -1;
		} else if (i == len - 1 && c != unknown && start != -1) {
			segmentList.push_back(make_pair(start, i));
			start = -1;
		}
	}
	// Post condition
	if (segmentList.empty()) {
		cerr << "KmerHistogram: At least one valid segment is required.";
		cerr << endl;
		cerr << "Sequence: " << *sequence << std::endl;
		throw std::exception();
	}

	// The hashed values, i.e. the values of the histograms.
	// The index is the 4ry representation of the key
	V *valueList = new V[maxTableSize];
	// Initialize values
	for (I i = 0; i < maxTableSize; i++) {
		valueList[i] = 0;
	}

	// Increment k-mer's in each valid segment
	for (auto segment : segmentList) {
		if (segment.first <= segment.second - k + 1) {
			vector < I > hashList = vector<I>();
			hashList.reserve(segment.second - segment.first + 1);
			hash(sequence, segment.first, segment.second - k + 1, &hashList);

			unsigned int size = hashList.size();
			for (unsigned int i = 0; i < size; i++) {
				valueList[hashList[i]]++;
			}
		}
	}

	// Check overflow
	for (I y = 0; y < maxTableSize; y++) {
		if (valueList[y] < 0) {
			cerr << "A negative value is a likely indication of overflow.";
			cerr << endl;
			cerr
					<< "To the developer: Consider larger data type in KmerHistogram.";
			cerr << endl;
			throw std::exception();
		}
	}

	return valueList;
}

/**
 * Make a list of the k-mers.
 */
template<class I, class V>
void KmerHistogram<I, V>::getKeys(vector<string> &keyList) {
	// Pre-condition: keyList must be empty.
	if (!keyList.empty()) {
		std::cerr << "Error at getKeys. ";
		std::cerr << "The output vector must be empty." << std::endl;
		throw std::exception();
	}
	keyList.reserve(maxTableSize);

	vector < string > alpha;
	vector < string > wordList;

	if (Parameters::isDNA()) {
		alpha = { "C", "T", "A", "G" };
		wordList = { "C", "T", "A", "G" };
	} else {
		std::cerr << "Error at getKeys. ";
		std::cerr << "Modes other than DNA are not supported yet.";
		std::cerr << std::endl;
		throw std::exception();
	}

	for (int i = 1; i < k; i++) {
		vector < string > temp;
		for (int j = 0; j < wordList.size(); j++) {
			for (int h = 0; h < alpha.size(); h++) {
				temp.push_back(wordList[j] + alpha[h]);
			}
		}
		wordList.clear();
		wordList = temp;
	}

	for (string s : wordList) {
		keyList.push_back(s);
	}
}

template<class I, class V>
void KmerHistogram<I, V>::getKeysDigitFormat(uint8_t keyList[]) {
	int alphaSize = Parameters::getAlphabetSize();

	for (int c = k - 1; c >= 0; c--) {
		for (int r = 0; r < maxTableSize; r++) {
			keyList[(r * k) + k - 1 - c] = (r / ((uint64_t) pow(alphaSize, c)))
					% alphaSize;
		}
	}
}

/**
 * Print the contents of the whole table
 */
template<class I, class V>
void KmerHistogram<I, V>::printTable(V *valueList) {
	vector < string > keys;
	getKeys (keys);
	for (I i = 0; i < keys.size(); i++) {
		string kmer = keys.at(i);
		cout << kmer << " -> " << (int) valueList[hash(&kmer)] << endl;
	}

	keys.clear();
}

/**
 * Print the contents of the histogram in Python's array format
 */
template<class I, class V>
void KmerHistogram<I, V>::printPythonFormat(V *valueList) {
	cout << "[";
	for (I i = 0; i < maxTableSize; i++) {
		cout << (uint64_t) valueList[i];
		if (i != maxTableSize - 1) {
			cout << ",";
		}
	}
	cout << "]" << endl;
}

template<class I, class V>
int KmerHistogram<I, V>::getK() {
	return k;
}

template<class I, class V>
I KmerHistogram<I, V>::getMaxTableSize() {
	return maxTableSize;
}
