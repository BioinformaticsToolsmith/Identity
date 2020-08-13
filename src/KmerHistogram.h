/*
 * KmerHistogram.h
 *
 *  Created on: Jul 25, 2012
 *  	Modified on: 12/21/2019
 *      Author: Hani Zakaria Girgis, PhD - NCBI/NLM/NIH
 */

#ifndef KMERHISTOGRAM_H_
#define KMERHISTOGRAM_H_

#include <string>
#include <vector>
// Review these include statement
#include <stdio.h>
#include <math.h>
#include <iostream>

#include "Parameters.h"

using namespace std;

template<class I, class V>
class KmerHistogram {

protected:
	/* Fields */
	static const int maxKeyLength = 15;
	const int k;
	I maxTableSize;

private:
	// [4^0, 4^1, ... , 4^(k-1)]
	I *bases;
	I *mMinusOne;
	void initialize(int, V);

	int digitList['T'+1];

public:
	/* Methods */
	KmerHistogram(int);
	virtual ~KmerHistogram();

	I hash(const string*);
	I hash(const string*, int);
	void hash(const string*, int, int, vector<I>*);
	V* build(const string *sequence);

	void getKeys(vector<string> &keys);
	void getKeysDigitFormat(uint8_t keyList[]);
	void printTable(V*);
	void printPythonFormat(V*);

	int getK();
	I getMaxTableSize();

};

#include "KmerHistogram.cpp"

#endif /* KMERHISTOGRAM_H_ */
