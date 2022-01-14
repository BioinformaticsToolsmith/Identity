
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
 * FeatureListSerializer.cpp
 *
 *  Created on: Jun 2, 2021
 *      Author: Hani Z. Girgis, PhD
 */

#include "Serializer.h"

std::ostream& operator<<(std::ostream &out, Feature &f) {
	out << std::setprecision(16);
	out << f.getNumOfComp() << "\t";
	if (f.getNumOfComp() == 0) {
		out << f.getFunIndex() << "\t";
	} else {
		out << -1 << "\t";
	}
	out << Serializer::changeSpaceToColumn(f.getName()) << "\t";
	out << f.getIsDistance() << "\t";
	out << f.getIsNormalized() << "\t";
	out << f.getNormP1() << "\t";
	out << f.getNormP2() << "\t";
	out << f.getTableIndex() << "\t";
	out << f.getIsSelected() << "\t";
	out << f.getIsNeeded() << "\t";
	out << f.getIsConverted() << "\t";
	out << f.getW() << std::endl;

	return out;
}

Serializer::Serializer(std::vector<Feature*> *featList, double *compList, int k,
		int histSize, double absError, int64_t maxLength, std::string file) {

	canOwnData = false;
	this->featList = featList;
	this->compList = compList;
	this->k = k;
	this->histSize = histSize;
	this->absError = absError;
	this->maxLength = maxLength;

	std::ofstream out(file);
	out << std::setprecision(16);
	out << k << std::endl;
	out << histSize << std::endl;
	out << absError << std::endl;
	out << maxLength << std::endl;

	for (int i = 0; i < Parameters::getAlphabetSize(); i++) {
		out << compList[i] << "\t";
	}
	out << std::endl;

	for (Feature *f : *featList) {
		out << (*f);
	}
	out.close();
}

Serializer::Serializer(std::string file) {
	canOwnData = true;
	std::ifstream in(file);

	// Fill k
	in >> k;
	// Fill histSize
	in >> histSize;
	// Fill absError;
	in >> absError;
	// Fill maxLength;
	in >> maxLength;

	// Fill composition list
	compList = new double[Parameters::getAlphabetSize()] { 0.0 };
	in >> compList[0] >> compList[1] >> compList[2] >> compList[3];

	// Fill feature list
	featList = new std::vector<Feature*>();
	int compNum;
	int funIndex;
	std::string name;
	bool isDistance;
	bool isNormalized;
	double normP1;
	double normP2;
	int tableIndex;
	bool isSelected;
	bool isNeeded;
	bool isConverted;
	double w;

	std::map<std::string, Feature*> nameFeatureMap;

	while (in >> compNum >> funIndex >> name >> isDistance >> isNormalized
			>> normP1 >> normP2 >> tableIndex >> isSelected >> isNeeded
			>> isConverted >> w) {
		// Make a new feature based on the number of components
		Feature *f;
		if (compNum == 0) {
			f = new Feature(funIndex, name, isDistance);
			nameFeatureMap[name] = f;
		} else if (compNum == 1) {
			std::string sName = Serializer::extractSquaredFeatureName(name);
			auto itr = nameFeatureMap.find(sName);
			if (itr != nameFeatureMap.end()) {
				f = new FeatureSquared(itr->second);
				nameFeatureMap[name] = f;
			} else {
				std::cerr << "There is no feature with the name of: " << sName
						<< std::endl;
				throw std::exception();
			}
		} else if (compNum == 2) {
			auto p = Serializer::extractPairedFeatureName(name);
			auto itr1 = nameFeatureMap.find(p.first);
			auto itr2 = nameFeatureMap.find(p.second);

			if (itr1 == nameFeatureMap.end()) {
				std::cerr << "There is no feature with the name of: " << p.first
						<< std::endl;
				throw std::exception();
			}

			if (itr2 == nameFeatureMap.end()) {
				std::cerr << "There is no feature with the name of: "
						<< p.second << std::endl;
				throw std::exception();
			}

			f = new FeaturePaired(itr1->second, itr2->second);
		} else {
			std::cerr << "Invalid number of components: " << compNum
					<< std::endl;
			throw std::exception();
		}

		f->setIsNormalized(isNormalized);
		f->setNormP1(normP1);
		f->setNormP2(normP2);
		f->setTableIndex(tableIndex);
		if (isSelected) {
			f->setIsSelected();
		}
		if (isNeeded) {
			f->setIsNeeded();
		}
		f->setIsConverted(isConverted);
		f->setW(w);
		featList->push_back(f);
	}

	in.close();
}

Serializer::~Serializer() {
	if (canOwnData) {
		delete[] compList;

		for (Feature *fPtr : *featList) {
			delete fPtr;
		}
		featList->clear();
		delete featList;
	}
}

double* Serializer::getCompList() const {
	return compList;
}

std::vector<Feature*>* Serializer::getFeatList() const {
	return featList;
}

int Serializer::getHistSize() const {
	return histSize;
}

int Serializer::getK() const {
	return k;
}

double Serializer::getAbsError() const {
	return absError;
}

int64_t Serializer::getMaxLength() const {
	return maxLength;
}

