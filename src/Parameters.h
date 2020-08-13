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
 * Parameters.h
 *
 *  Created on: Nov 30, 2019
 *      Author: Hani Z. Girgis, PhD
 *      Purpose: Holds all parameter values as static variables, so they can be accessed any where in the program
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <iostream>
#include <limits>

class Parameters {
private:
	static int MODE;

	// Parameters controlling the trainer.
	// static double THRESHOLD;
	static double MIN_ID;
	static int MUT_PER_TEMP;
	static int BLOCK_SIZE;

	// Parameters controlling sequence generation
	static bool IS_SINGLE;
	static bool IS_BLOCK;
	static bool IS_TRANSLOCATION;
	static bool IS_INVERSTION;
	static int MIN_BLOCK_SIZE;
	static int MAX_BLOCK_SIZE;

	// Parameters controlling training the classifier and the regression model
	static int MIN_FEAT_NUM;
	static int MAX_FEAT_NUM;

	static double DELTA_R;
	static double DELTA_C;

	static int K_RELAX;
	static void checkMutPerTemp();

public:
	static const int DNA = 0;
	static const int PROTEIN = 1;
	Parameters();
	virtual ~Parameters();
	static int getMode();
	static void setMode(int);
	static char getUnknown();
	static int getAlphabetSize();
	static bool isDNA();

	// Parameters controlling the classifier
	// static double getThreshold();
	// static void setThreshold(double);
	// static void checkThreshold();

	static double getMinId();
	static void setMinId(double);
	static void checkMinId();

	static int getMutPerTemp();
	static void setMutPerTemp(int);

	// Parameters controlling the trainer
	static int getBlockSize();
	static void setBlockSize(int);

	// Parameters controlling sequence generation
	static void enableMutSingle();
	static void enableMutBlock();
	static void enableMutInverstion();
	static void enableMutTranslocation();
	static bool getMutSingle();
	static bool getMutBlock();
	static bool getMutInverstion();
	static bool getMutTranslocation();

	static int getMutMinBlockSize();
	static int getMutMaxBlockSize();
	static void setMutMinBlockSize(int);
	static void setMutMaxBlockSize(int);

	static int getKRelax();
	void setKRelax(int);

	// Parameters controlling training
	static int getMinFeatNum();
	static void setMinFeatNum(int);
	static int getMaxFeatNum();
	static void setMaxFeatNum(int);

	static double getDeltaC();
	static void setDeltaC(double deltaC);
	static double getDeltaR();
	static void setDeltaR(double deltaR);
};

#endif /* PARAMETERS_H_ */
