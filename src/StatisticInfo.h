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
 * StatisticInfo.h
 *
 *  Created on: May 19, 2020
 *      Author: Dr. Hani Zakaria Girgis
 */

#ifndef STATISTICINFO_H_
#define STATISTICINFO_H_

#include<vector>
#include "Statistician.h"
#include "Feature.h"

/**
 * The single instance of this class keeps information
 * about the statistics implemented by the statistician
 */
class StatisticInfo {
private:
	std::vector<Feature*> *fList;
	static StatisticInfo *singleton;

	StatisticInfo();

public:
	static StatisticInfo* getInstance();
	virtual ~StatisticInfo();
	std::vector<Feature*>* getList();
	int getStatNum();
};

#endif /* STATISTICINFO_H_ */
