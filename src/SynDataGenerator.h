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
 * Trainer.h
 *
 *  Created on: Feb 25, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef SYNDATAGENERATOR_H_
#define SYNDATAGENERATOR_H_

#include "KmerHistogram.h"
#include "Matrix.h"
#include "Statistician.h"
#include "Mutator.h"
#include "DataGenerator.h"
#include "StatisticInfo.h"

class SynDataGenerator: public DataGenerator {
private:
	double threshold;
	std::vector<int> funIndexList;
	int threadNum;

	template<class V>
	void generateDataHelper();
	void fillCompositionList();

public:
	SynDataGenerator(std::string, double, int, std::vector<int> funIndexList =
			std::vector<int>());
	SynDataGenerator(std::string, std::string, double, int,
			std::vector<int> funIndexList = std::vector<int>());
	virtual ~SynDataGenerator();
	virtual void generateData();
	virtual void clearData();
};

#endif /* TRAINERTEMPLATE_H_ */
