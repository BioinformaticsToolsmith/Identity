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
 * FastaReader.h
 *
 *  Created on: Nov 2, 2019
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef FASTAREADER_H_
#define FASTAREADER_H_

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm> // for_each
#include <stdio.h> // fread

#include "Parameters.h"

using namespace std;

typedef vector<pair<string*, string*> > Block;

class FastaReader {
private:
	int blockSize;
	string fileName;
	ifstream in;
	int maxLen;

	bool isDone;
	long int currentPos;
	char unknown;
	char codeMap[128];
	const char NOT = '!';

	bool isAllInvalid(const string*);

public:
	FastaReader(string, int, long int currentPosIn = 0, int maxLenIn = 0);
	virtual ~FastaReader();
	Block* read();
	//Block * fRead();
	static void deleteBlock(Block*);
	bool isStillReading();
	long int getCurrentPos();
	int getMaxLen();
	void restart();
	void setBlockSize(int);
};

#endif /* FASTAREADER_H_ */
