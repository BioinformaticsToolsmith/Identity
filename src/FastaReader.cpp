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
 * FastaReader.cpp
 *
 *  Created on: Nov 2, 2019
 *      Author: Hani Zakaria Girgis, PhD
 *
 *
 *  Notes:
 *  + A sequence of all unknown nucleotides or aa's is excluded.
 *  + A sequence is converted to upper case.
 */

#include "FastaReader.h"

FastaReader::FastaReader(std::string fileNameIn, int blockSizeIn,
		long int currentPosIn, int maxLenIn) {
	fileName = fileNameIn;
	blockSize = blockSizeIn;

	unknown = Parameters::getUnknown();

	isDone = false;
	currentPos = currentPosIn;
	in = std::ifstream(fileName.c_str());
	if (!in.good()) {
		std::cerr << "Cannot open file: " << fileName << std::endl;
		throw std::exception();
	}

	in.seekg(currentPos);
	if (in.rdstate() & (in.badbit | in.eofbit | in.failbit)) {
		std::cerr << "FastaReader Error: " << std::endl;
		std::cerr << "Cannot move to position: " << currentPos << std::endl;
		throw std::exception();
	}

	maxLen = maxLenIn;

	for (int i = 0; i < 128; i++) {
		codeMap[i] = NOT;
	}

	codeMap['A'] = 'A';
	codeMap['C'] = 'C';
	codeMap['G'] = 'G';
	codeMap['T'] = 'T';
	codeMap['N'] = 'N';
	codeMap['R'] = 'G';
	codeMap['Y'] = 'C';
	codeMap['M'] = 'A';
	codeMap['K'] = 'T';
	codeMap['S'] = 'G';
	codeMap['W'] = 'T';
	codeMap['H'] = 'C';
	codeMap['B'] = 'T';
	codeMap['V'] = 'A';
	codeMap['D'] = 'T';
	// Added on 11/19/2020 to enable reading multi-sequence alignments
	codeMap['-'] = '-';
}

FastaReader::~FastaReader() {
	in.close();
}

/**
 * Utility function
 */
void FastaReader::deleteBlock(Block *block) {
	for (auto p : *block) {
		delete p.first;
		delete p.second;
	}

	block->clear();
	// Added on 4/25/2020.
	delete block;
}

Block* FastaReader::read() {
	bool isFirst = true;
	string *base = new string("");
	string *info;
	Block *b = new Block();
	b->reserve(blockSize);

	int counter = 0;
	while (in.good() && counter < blockSize) {
		string line;
		getline(in, line);
		int len = line.length();

		if(line[len-1] == '\r'){
			len--;
			line.pop_back();
		}

		if (line[0] == '>') {
			if (!isFirst) {
				if (base->length() > maxLen) {
					maxLen = base->length();
				}

				if (isAllInvalid(base)) {
					delete info;
					delete base;
				} else {
					base->shrink_to_fit();
					b->push_back(make_pair(info, base));
					counter++;
				}

				if (counter < blockSize) {
					info = new string(line);
					base = new string("");
					base->reserve(maxLen);
				} else {
					// Go back one line
					in.seekg(-1 - len, in.cur);
				}
			} else {
				info = new string(line);
				isFirst = false;
			}
		} else {
			// Convert non-traditional bases to traditional ones
			for (int h = 0; h < len; h++) {
				// Convert a char to upper case if needed
				char & o = line[h];
				if (o >= 97) {
					line[h] = o - 32;
				}

				char c = codeMap[o];
				if (c != NOT) {
					line[h] = c;
				} else {
					std::cerr << "Something wrong with: " << *info << std::endl;
					std::cerr << "At this line: " << line << std::endl;
					std::cerr << "Invalid nucleotide symbol: ("
							<< (char) line[h];
					std::cerr << ")" << endl;
					throw std::exception();
				}
			}

			base->append(line);
		}
	}
	// Record current position
	currentPos = in.tellg();

	// The file ended and the last sequence has not been added to
	// the block
	if (!in.good()) {
		if (isAllInvalid(base)) {
			delete info;
			delete base;
		} else {
			base->shrink_to_fit();
			b->push_back(make_pair(info, base));
		}
		isDone = true;
		b->shrink_to_fit();
	}

	return b;
}

/**
 * Check if the sequence consists of all unknown nucleotides or aa's
 */
bool FastaReader::isAllInvalid(const string *base) {
	bool r = true;
	const char *arr = base->c_str();
	int len = base->length();
	for (int i = 0; i < len; i++) {
		if (arr[i] != unknown) {
			r = false;
			break;
		}
	}

	return r;
}

bool FastaReader::isStillReading() {
	return !isDone;
}

long int FastaReader::getCurrentPos() {
	return currentPos;
}

int FastaReader::getMaxLen() {
	return maxLen;
}
