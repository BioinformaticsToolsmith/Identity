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
 * StatisticInfo.cpp
 *
 *  Created on: May 19, 2020
 *      Author: Dr. Hani Zakaria Girgis
 */

#include "StatisticInfo.h"

StatisticInfo *StatisticInfo::singleton = nullptr;

StatisticInfo::StatisticInfo() {

	fList = new std::vector<Feature*>();
	/**
	 * These are distance statistics
	 */
	fList->push_back(new Feature(Stat::MANHATTAN, "manhattan", true));
	fList->push_back(new Feature(Stat::EUCLIDEAN, "euclidean", true));
	fList->push_back(new Feature(Stat::CHI_SQUARED, "chi_squared", true));
	fList->push_back(new Feature(Stat::CHEBYSHEV, "chebyshev", true));
	fList->push_back(new Feature(Stat::HAMMING, "hamming", true));
	fList->push_back(new Feature(Stat::MINKOWSKI, "minkowski", true));
	fList->push_back(new Feature(Stat::COSINE, "cosine", true));
	fList->push_back(new Feature(Stat::CORRELATION, "correlation", true));
	fList->push_back(new Feature(Stat::BRAYCURTIS, "bray_curtis", true));
	fList->push_back(new Feature(Stat::SQUARED_CHORD, "squared_chord", true));
	fList->push_back(new Feature(Stat::HELLINGER, "hellinger", true));
	fList->push_back(
			new Feature(Stat::CUMULATIVE_DIFF, "cumulative_difference", true));
	fList->push_back(new Feature(Stat::EMD, "emd", true));
	fList->push_back(new Feature(Stat::KL_CONDITIONAL, "kl_conditional", true));
	fList->push_back(new Feature(Stat::K_DIVERGENCE, "k_divergence", true));
	fList->push_back(
			new Feature(Stat::JEFFREY_DIVERGENCE, "jeffrey_divergence", true));
	fList->push_back(
			new Feature(Stat::JENSEN_SHANNON_DIVERGENCE,
					"jensen_shannon_divergence", true));
	fList->push_back(new Feature(Stat::RRE, "rre", true));

	/**
	 * These are similarity statistics
	 */
	fList->push_back(new Feature(Stat::INTERSECTION, "intersection", false));
	fList->push_back(new Feature(Stat::KULCZYNSKI_1, "kulczynski_1", false));
	fList->push_back(new Feature(Stat::KULCZYNSKI_2, "kulczynski_2", false));
	fList->push_back(new Feature(Stat::COVARIANCE_R, "covariance_r", false));
	fList->push_back(
			new Feature(Stat::HARMONIC_MEAN_R, "harmonic_mean_r", false));
	fList->push_back(new Feature(Stat::SIM_RATIO, "sim_ratio", false));
	fList->push_back(new Feature(Stat::MARKOV_R, "markov_r", false));
	fList->push_back(new Feature(Stat::SIM_MM, "simMM", false));
	fList->push_back(new Feature(Stat::LENGTH_RATIO, "length_ratio", false));
	fList->push_back(new Feature(Stat::D2S_R, "d2_s_r", false));
	fList->push_back(new Feature(Stat::D2STAR, "d2_star", false));
}

StatisticInfo::~StatisticInfo() {
	for (auto f : *fList) {
		delete f;
	}
	fList->clear();
	delete fList;
}

StatisticInfo* StatisticInfo::getInstance() {
	if (singleton == nullptr) {
		singleton = new StatisticInfo();
	}
	return singleton;
}

std::vector<Feature*>* StatisticInfo::getList() {
	return fList;
}

int StatisticInfo::getStatNum() {
	return fList->size();
}
