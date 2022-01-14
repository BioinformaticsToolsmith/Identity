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
 * FeatureListSerializer.h
 *
 *  Created on: Jun 2, 2021
 *      Author: Dr. Hani Z. Girgis
 */

#ifndef SRC_SERIALIZER_H_
#define SRC_SERIALIZER_H_

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <map>

#include "Feature.h"
#include "FeatureSquared.h"
#include "FeaturePaired.h"
#include "Parameters.h"

class Serializer {
private:
	std::vector<Feature*> *featList;
	double *compList;
	int k;
	int histSize;
	double absError;
	bool canOwnData;
	int64_t maxLength;

public:
	Serializer(std::vector<Feature*>*, double*, int, int, double, int64_t,
			std::string);
	Serializer(std::string);

	virtual ~Serializer();

	static inline std::string changeSpaceToColumn(std::string in) {
		for (int i = 0; i < in.size(); i++) {
			if (in.at(i) == ' ') {
				in.at(i) = ':';
			}
		}
		return in;
	}

	// kulczynski_2^2
	static inline std::string extractSquaredFeatureName(std::string in) {
		std::string name("");
		for (int i = 0; i < in.size(); i++) {
			if (in[i] == '^') {
				break;
			} else {
				name.append(1, in[i]);
			}
		}
		return name;
	}

	// euclidean:x:sim_ratio
	static inline std::pair<std::string, std::string> extractPairedFeatureName(
			std::string in) {
		std::string first("");
		int i;
		for (i = 0; i < in.size(); i++) {
			if (in[i] == ':') {
				break;
			} else {
				first.append(1, in.at(i));
			}
		}
		i++; // x
		i++; // :
		std::string second("");
		for (i++; i < in.size(); i++) {
			second.append(1, in.at(i));
		}
		return std::make_pair(first, second);
	}

	// void writeList(std::vector<Feature*> &featList, std::string file);
	// std::vector<Feature*>* readList(std::string file);

	double* getCompList() const;
	std::vector<Feature*>* getFeatList() const;
	int getHistSize() const;
	int getK() const;
	double getAbsError() const;
	int64_t getMaxLength() const;
};

#endif /* SRC_SERIALIZER_H_ */
