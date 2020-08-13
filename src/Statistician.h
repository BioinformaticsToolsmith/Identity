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
 * Statistician.h
 *
 *  Created on: Dec 21, 2019
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef STATISTICIAN_H_
#define STATISTICIAN_H_

#include <cmath> // pow, sqrt, round
#include <iostream> // cout
#include <algorithm>
#include <any>
#include <functional>
#include <fstream>
#include <vector>
#include <limits>

#include "Parameters.h"
#include "Util.h"
#include "Feature.h"

// Enumerator of all statistics
enum Stat : int8_t {
	/*
	 * Add new distance statistics after this block comment.
	 */
	MANHATTAN,
	EUCLIDEAN,
	CHI_SQUARED,
	CHEBYSHEV,
	HAMMING,
	MINKOWSKI,
	COSINE,
	CORRELATION,
	BRAYCURTIS,
	SQUARED_CHORD,
	HELLINGER,
	CUMULATIVE_DIFF,
	EMD,
	KL_CONDITIONAL,
	K_DIVERGENCE,
	JEFFREY_DIVERGENCE,
	JENSEN_SHANNON_DIVERGENCE,
	RRE,

	/*
	 * Add new distance statistics before this block comment.
	 */
	Dist_NUM,

	/**
	 * Add new similarity statistics after this block comment.
	 */
	INTERSECTION,
	KULCZYNSKI_1,
	KULCZYNSKI_2,
	COVARIANCE_R,
	HARMONIC_MEAN_R,
	SIM_RATIO,
	MARKOV_R,
	SIM_MM,
	LENGTH_RATIO,
	D2S_R,
	D2STAR,

	/**
	 * Add new similarity statistics before this block comment.
	 */
	ALL_NUM
};

// Important: V must be signed integer type!
template<class V>
class Statistician {
private:
	const int histogramSize;
	const int k;
	const V *h1; // Kmer histogram of sequence 1
	const V *h2; // Kmer histogram of sequence 2
	const uint64_t *mono1; // Monomer histogram of sequence 1
	const uint64_t *mono2; // Monomer histogram of sequence 2
	// Array representing a background model for C, T, A and ,G, e.g. n[4] = {0.25, 0.25, 0.25, 0.25}.
	const double *background;
	const uint8_t *keyList; // key list in digit format, e.g. getKeysDigitFormat defined in KmerHistogram

	double mean1; // Mean of kmer histogram 1
	double mean2; // Mean of kmer histogram 2
	double *p1; // Probability vector based on kmer histogram 1
	double *p2; // Probability vector based on kmer histogram 2
	V *mean1And2;

	// methodList is an array of function pointers
	static double (Statistician<V>::*methodList[Stat::ALL_NUM])();

public:
	Statistician(int histogramSizeIn, int kIn, const V *h1In, const V *h2In,
			const uint64_t *mono1In, const uint64_t *mono2In,
			const double *backgroundIn, const uint8_t*);
	virtual ~Statistician();

	/**
	 * Distance statistics
	 */

	double manhattanDistance();

	double euclideanDistance();

	double canberraDistance();

	double chiSquaredDistance();

	double chebyshevDistance();

	double hammingDistance();

	double minkowskiDistance();

	double cosineDistance();

	double correlationDistance();

	double braycurtisDistance();

	double squaredChordDistance();

	double hellingerDistance();

	double cumulativeDiffDistance();

	double emdDistance();

	double afdDistance();

	double klConditionalDistance();

	double kDivergenceDistance();

	double jeffreyDivergenceDistance();

	double klDivergenceSymmetricDistance();

	double jensenShannonDivergenceDistance();

	double jensenShannonDivergenceGDistance();

	double jensenShannonDivergenceHDistance();

	double rreDistance();

	/**
	 * Similarity statistics
	 */

	double intersectionSimilarity();

	double kulczynski1Similarity();

	double kulczynski2Similarity();

	double covarianceSimilarity();

	double covarianceRSimilarity();

	double harmonicMeanSimilarity();

	double harmonicMeanRSimilarity();

	double simRatioSimilarity();

	double markovSimilarity();

	double markovRSimilarity();

	double simMMSimilarity();

	double lengthRatioSimilarity();

	double d2sSimilarity();

	double d2sRSimilarity();

	double d2starSimilarity();

	/**
	 * Helper methods
	 */
	double cosineDistanceHelper(const V*, const V*);
	double klDivergenceDistanceHelper(const double*, const double*);
	double covarianceSimilarityHelper(const V*, const V*, double, double);
	double harmonicMeanSimilarityHelper(const V*, const V*);
	double markovSimilarityHelper(const V*, const V*);
	double d2sSimilarityHelper(const V*, const V*);

	/**
	 * Requires a vector of standard deviations calculated on the entire data set.
	 * Cannot be saved at this time in the function array.
	 */
	double seuclideanDistance(const double*);

	/**
	 * Calculate all statistics or subset of them
	 */
	void calculateAll(std::vector<double>& /*out*/);
	void calculate(std::vector<int>& /*in*/, std::vector<double>& /*out*/);
	void calculate(const int *s/*method index list*/,
			int size /*size of method index list*/, double *r /*result list*/);

	/**
	 * Inline methods/functions
	 */
	inline double mean(const V *h) {
		double m = 0.0;
		for (int i = 0; i < histogramSize; i++) {
			m += h[i];
		}
		return m / histogramSize;
	}

	inline uint64_t sum(const V *h) {
		uint64_t s = 0;
		for (int i = 0; i < histogramSize; i++) {
			s += h[i];
		}
		return s;
	}

	inline uint64_t sum(const uint64_t *h, int size) {
		uint64_t s = 0;
		for (int i = 0; i < size; i++) {
			s += h[i];
		}
		return s;
	}

	inline double norm(const V *h) {
		uint64_t n = 0;
		for (int i = 0; i < histogramSize; i++) {
			n += h[i] * h[i];
		}

		if (n < 0.0) {
			std::cerr << "Vector norm cannot be negative." << std::endl;
			throw std::exception();
		}

		return sqrt(n);
	}

	inline V absolute(V n) {
		return n < 0 ? -1.0 * n : n;
	}
};

#include "Statistician.cpp"

#endif /* STATISTICIAN_H_ */
