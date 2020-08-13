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
 * Feature.h
 *
 *  Created on: May 19, 2020
 *      Author: Hani Z. Girgis, PhD
 */

#ifndef FEATURE_H_
#define FEATURE_H_

#include <string>
#include <iostream>

class Feature {
private:
	// Function index that is used by the statistician to calculate it.
	const int funIndex;

	int tableIndex;

	// Distance or similarity
	const bool isDistance;
	bool isConverted;
	double w;

protected:
	// Feature name
	const std::string name;

	// These three parameters are needed to normalize this feature
	bool isNormalized;
	double normP1;
	double normP2;

	// These two parameters are needed by the feature selector
	bool isNeeded;
	bool isSelected;

public:
	Feature(int, std::string, bool);
	virtual ~Feature();
	Feature(const Feature &other);
	Feature(Feature &&other);

	// Methods for normalizing/standardizing the feature
	bool getIsNormalized() const;
	void setIsNormalized(bool isNormalized);
	double getNormP1() const;
	void setNormP1(double normP1);
	double getNormP2() const;
	void setNormP2(double normP2);

	int getTableIndex() const;
	void setTableIndex(int tableIndex);

	const std::string& getName() const;

	virtual const int getFunIndex() const;

	// Methods for composition
	virtual int getNumOfComp();
	virtual int getCompOneIndex();
	virtual int getCompTwoIndex();

	bool getIsDistance() const;

	bool getIsNeeded() const;
	virtual void setIsNeeded();
	bool getIsSelected() const;
	virtual void setIsSelected();

	virtual void setCompOne(Feature*);
	virtual void setCompTwo(Feature*);
	bool getIsConverted() const;
	void setIsConverted(bool isConverted);
	double getW() const;
	void setW(double w);
};

#endif /* FEATURE_H_ */
