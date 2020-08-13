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
 * Statistician.cpp
 *
 *  Created on: Dec 21, 2019
 *      Author: Hani Zakaria Girgis, PhD
 *
 */

// ToDo: Fix the issue with AFD
template<class V>
double (Statistician<V>::*Statistician<V>::methodList[Stat::ALL_NUM])() = {
	&Statistician < V > ::manhattanDistance,
	&Statistician < V > ::euclideanDistance,
	&Statistician < V> ::chiSquaredDistance,
	&Statistician < V > ::chebyshevDistance,
	&Statistician < V > ::hammingDistance,
	&Statistician < V > ::minkowskiDistance,
	&Statistician < V > ::cosineDistance,
	&Statistician < V > ::correlationDistance,
	&Statistician < V > ::braycurtisDistance,
	&Statistician < V > ::squaredChordDistance,
	&Statistician < V > ::hellingerDistance,
	&Statistician < V> ::cumulativeDiffDistance,
	&Statistician < V > ::emdDistance,
	&Statistician < V> ::klConditionalDistance,
	&Statistician < V > ::kDivergenceDistance,
	&Statistician < V > ::jeffreyDivergenceDistance,
	&Statistician < V > ::jensenShannonDivergenceDistance,
	&Statistician < V > ::rreDistance,
	nullptr,
	&Statistician < V > ::intersectionSimilarity,
	&Statistician < V > ::kulczynski1Similarity,
	&Statistician < V > ::kulczynski2Similarity,
	&Statistician < V > ::covarianceRSimilarity,
	&Statistician < V > ::harmonicMeanRSimilarity,
	&Statistician < V > ::simRatioSimilarity,
	&Statistician < V > ::markovRSimilarity,
	&Statistician < V > ::simMMSimilarity,
	&Statistician < V > ::lengthRatioSimilarity,
	&Statistician < V > ::d2sRSimilarity,
	&Statistician < V > ::d2starSimilarity
};

template<class V>
Statistician<V>::Statistician(int histogramSizeIn, int kIn, const V *h1In,
		const V *h2In, const uint64_t *mono1In, const uint64_t *mono2In,
		const double *backgroundIn, const uint8_t *keyListIn) :
		histogramSize(histogramSizeIn), k(kIn), h1(h1In), h2(h2In), mono1(
				mono1In), mono2(mono2In), background(backgroundIn), keyList(
				keyListIn) {

	// Calculate means
	mean1 = mean(h1);
	mean2 = mean(h2);

	if (Util::isEqual(mean1, 0.0) || Util::isEqual(mean2, 0.0)) {
		std::cerr << "Mean 1 (mean1) and Mean 2 (mean2) cannot be zeros.";
		std::cerr << std::endl;
		throw std::exception();
	}

	// Calculate probability vectors with pseudo counts.
	uint64_t s1 = sum(h1) + histogramSize;
	uint64_t s2 = sum(h2) + histogramSize;
	p1 = new double[histogramSize];
	p2 = new double[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		p1[i] = (h1[i] + 1.0) / s1;
		p2[i] = (h2[i] + 1.0) / s2;
	}

	// Calculate mean vector element wise
	mean1And2 = new V[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		uint64_t m = h1[i] + h2[i];
		mean1And2[i] = round(m / 2.0);
	}
}

template<class V>
Statistician<V>::~Statistician() {
	delete[] p1;
	delete[] p2;
	delete[] mean1And2;
}

template<class V>
double Statistician<V>::manhattanDistance() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d += absolute(h1[i] - h2[i]);
	}
	return d;
}

template<class V>
double Statistician<V>::euclideanDistance() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		double temp = h1[i] - h2[i];
		d += temp * temp;
	}
	return sqrt(d);
}

template<class V>
double Statistician<V>::canberraDistance() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		// Skip row if both entries are zeros.
		if (h1[i] > 0 || h2[i] > 0) {
			d += absolute(h1[i] - h2[i]) / ((double) h1[i] + h2[i]);
		}
	}
	return d;
}

template<class V>
double Statistician<V>::chiSquaredDistance() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		// Skip row if both entries are zeros.
		if (h1[i] > 0 || h2[i] > 0) {
			double diff = h1[i] - h2[i];
			d += (diff * diff) / (h1[i] + h2[i]);
		}
	}
	return d;
}

template<class V>
double Statistician<V>::chebyshevDistance() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		V diff = absolute(h1[i] - h2[i]);
		if (diff > d) {
			d = diff;
		}
	}
	return d;
}

template<class V>
double Statistician<V>::hammingDistance() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		if (h1[i] != h2[i]) {
			d++;
		}
	}
	return d / histogramSize;
}

/**
 * Works only on power of 3
 */
template<class V>
double Statistician<V>::minkowskiDistance() {
	long long int d = 0.0;

	for (int i = 0; i < histogramSize; i++) {
		V z = absolute(h1[i] - h2[i]);
		d += (z * z * z);
	}

	return std::cbrt(d);
}

/**
 * Helper function to cosineDistance() and correlationDistance()
 */
template<class V>
double Statistician<V>::cosineDistanceHelper(const V *v1, const V *v2) {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d += v1[i] * v2[i];
	}

	double n1 = norm(v1);
	double n2 = norm(v2);
	if (Util::isEqual(n1, 0.0) || Util::isEqual(n2, 0.0)) {
		std::cerr << "Error at Cosine Distance. ";
		std::cerr << "Vector norm (n1 or n2) is zero." << std::endl;
		throw std::exception();
	}

	return 1.0 - d / (n1 * n2);
}

template<class V>
double Statistician<V>::cosineDistance() {
	return cosineDistanceHelper(h1, h2);
}

template<class V>
double Statistician<V>::correlationDistance() {
	// Center h1 around its mean.
	V *n1 = new V[histogramSize];
	double m1 = round(mean1);
	for (int i = 0; i < histogramSize; i++) {
		n1[i] = h1[i] - m1;
	}

	// Center h2 around its mean.
	V *n2 = new V[histogramSize];
	double m2 = round(mean2);
	for (int i = 0; i < histogramSize; i++) {
		n2[i] = h2[i] - m2;
	}

	// Calculate the cosine distance on the centered histograms.
	double r = cosineDistanceHelper(n1, n2);

	delete[] n1;
	delete[] n2;

	return r;
}

template<class V>
double Statistician<V>::braycurtisDistance() {
	double d1 = 0.0;
	double d2 = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d1 += absolute(h1[i] - h2[i]);
		d2 += h1[i] + h2[i];
	}

	if (Util::isEqual(d2, 0.0)) {
		std::cerr << "Error at Bray-curtis distance. ";
		std::cerr << "The denominator (d2) is zero";
		throw std::exception();
	}

	return d1 / d2;
}

template<class V>
double Statistician<V>::squaredChordDistance() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d += h1[i] + h2[i] - 2 * sqrt(h1[i] * h2[i]);
	}
	return d;
}

template<class V>
double Statistician<V>::hellingerDistance() {
	// Divide h1 by its mean.
	double *n1 = new double[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		n1[i] = h1[i] / mean1;
	}

	// Divide h2 by its mean.
	double *n2 = new double[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		n2[i] = h2[i] / mean2;
	}

	// Calculate Squared Chord Distance on n1 & n2
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d += n1[i] + n2[i] - 2 * sqrt(n1[i] * n2[i]);
	}

	delete[] n1;
	delete[] n2;

	return sqrt(2 * d);
}

/**
 * Reference: "A Closed-form Gradient for the 1D Earth Mover’s Distance for
 * Spectral Deep Learning on Biological Data."
 * http://www.cs.toronto.edu/~makarand/papers/ICML2016_CompBioWorkshop.pdf
 */
template<class V>
double Statistician<V>::emdDistance() {
	double cumulativeDiff = 0.0;
	double emd = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		cumulativeDiff += p1[i] - p2[i];
		emd += std::abs(cumulativeDiff);
	}
	return emd;
}

/**
 * Based on code by Ben James.
 */
template<class V>
double Statistician<V>::cumulativeDiffDistance() {
	double c1 = 0.0, c2 = 0.0;
	double cumulativeDiff = 0.0;
	for (auto i = 0; i < histogramSize; i++) {
		c1 += h1[i];
		c2 += h2[i];
		cumulativeDiff += c1 > c2 ? c1 - c2 : c2 - c1;
	}
	return cumulativeDiff;
}

template<class V>
double Statistician<V>::afdDistance() {
	if (k != 2) {
		std::cerr << "AFD cannot be calculated for k other than 2: Received: ";
		std::cerr << k << std::endl;
		throw std::exception();
	}
	const auto alphaSize = Parameters::getAlphabetSize();

	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		double r = (double) h1[i] / (mono1[i % alphaSize] + 1.0)
				- (double) h2[i] / (mono2[i % alphaSize] + 1.0);
		double fm = 1.0 / pow(1.0 + r, -14);
		d += r * fm * r * fm;

		if (std::isinf(d)) {
			std::cerr << (double) h1[i] / (mono1[i % alphaSize] + 1.0) << " ";
			std::cerr << (double) h2[i] / (mono2[i % alphaSize] + 1.0) << " ";
			std::cerr << fm << " " << r * fm * r * fm << std::endl;
			throw std::exception();
		}
	}

	return d;
}

template<class V>
double Statistician<V>::kDivergenceDistance() {
	// This is an asymmetric statistics. So average two
	// statistics, each of which is with respect to one histogram.
	double d1 = 0.0;
	double d2 = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		double avg = (p1[i] + p2[i]) / 2.0;
		d1 += p1[i] * log(p1[i] / avg);
		d2 += p2[i] * log(p2[i] / avg);
	}
	return (d1 + d2) / 2.0;
}

template<class V>
double Statistician<V>::jeffreyDivergenceDistance() {
	// This is a symmetric statistic.
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d += (p1[i] - p2[i]) * log(p1[i] / p2[i]);
	}
	return d;
}

template<class V>
double Statistician<V>::klDivergenceDistanceHelper(const double *o1,
		const double *o2) {
	// This is an asymmetric statistic.
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d += o1[i] * log(o1[i] / o2[i]);
	}
	return d;
}

/**
 * Same as klDivergenceSymmeticDistance but works on probability vectors
 * not histograms.
 */
template<class V>
double Statistician<V>::klDivergenceSymmetricDistance() {
	return (klDivergenceDistanceHelper(p1, p2)
			+ klDivergenceDistanceHelper(p2, p1)) / 2.0;
}

template<class V>
double Statistician<V>::jensenShannonDivergenceDistance() {
	// The sum includes pseudo counts.
	uint64_t s = sum(mean1And2) + histogramSize;

	// Add pseudo count and convert to probability.
	double *p = new double[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		p[i] = (mean1And2[i] + 1.0) / s;
	}

	double r = (klDivergenceDistanceHelper(p1, p)
			+ klDivergenceDistanceHelper(p2, p)) / 2.0;

	delete[] p;

	return r;
}

template<class V>
double Statistician<V>::jensenShannonDivergenceGDistance() {
	// Array a is the geometric mean of h1 and h2.
	V *a = new V[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		uint64_t m = h1[i] * h2[i];
		a[i] = round(sqrt(m));
	}

	// The sum includes pseudo counts.
	uint64_t s = sum(a) + histogramSize;

	// Add pseudo count and convert to probability.
	double *p = new double[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		p[i] = (a[i] + 1.0) / s;
	}

	double r = (klDivergenceDistanceHelper(p1, p)
			+ klDivergenceDistanceHelper(p2, p)) / 2.0;

	delete[] a;
	delete[] p;

	return r;
}

template<class V>
double Statistician<V>::jensenShannonDivergenceHDistance() {
	// Array a is the harmonic mean of h1 and h2.
	V *a = new V[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		double m = 1.0 / h1[i] + 1.0 / h2[i];
		a[i] = round(2.0 / m);
	}

	// The sum includes pseudo counts.
	uint64_t s = sum(a) + histogramSize;

	// Add pseudo count and convert to probability.
	double *p = new double[histogramSize];
	for (int i = 0; i < histogramSize; i++) {
		p[i] = (a[i] + 1.0) / s;
	}

	double r = (klDivergenceDistanceHelper(p1, p)
			+ klDivergenceDistanceHelper(p2, p)) / 2.0;

	delete[] a;
	delete[] p;

	return r;
}

/**
 * Based on code by Ben James and Hani Z. Girgis
 */
template<class V>
double Statistician<V>::klConditionalDistance() {
	const int a = Parameters::getAlphabetSize();
	uint64_t sum4_1 = a, sum4_2 = a; // Sum for every 4 nucleotides or 22 a.a.
	double outer_sum_1 = 0.0, outer_sum_2 = 0.0; // Prior K-mer sum

	for (auto i = 0; i < histogramSize; i++) {
		sum4_1 += h1[i];
		sum4_2 += h2[i];

		if (i % a == a - 1) { //finished counting word, now compute probabilities
			double inner_sum_1 = 0.0;    // Sum of p(X|Y) * log(p(X|Y) / q(X|Y))
			double inner_sum_2 = 0.0;    // Sum of q(X|Y) * log(q(X|Y) / p(X|Y))
			for (int j = i - (a - 1); j <= i; j++) {
				double conditional_1 = (double) (h1[j] + 1.0) / sum4_1;
				double conditional_2 = (double) (h2[j] + 1.0) / sum4_2;
				double lg = log(conditional_1 / conditional_2);
				inner_sum_1 += conditional_1 * lg;
				inner_sum_2 += -1 * conditional_2 * lg;
			}
			outer_sum_1 += sum4_1 * inner_sum_1;
			outer_sum_2 += sum4_2 * inner_sum_2;
			// Reinitialize to pseudo count
			sum4_1 = a;
			sum4_2 = a;
		}
	}
	double left = outer_sum_1 / (sum(h1) + histogramSize);
	double right = outer_sum_2 / (sum(h2) + histogramSize);
	return (left + right) / 2.0;
}

/**
 * Should work on DNA and Proteins
 */
template<class V>
double Statistician<V>::rreDistance() {
	double d1 = 0.0, d2 = 0.0;
	const int alphaSize = Parameters::getAlphabetSize();

	for (auto i = 0; i < histogramSize; i += alphaSize) {
		double sum1 = alphaSize, sum2 = alphaSize; // Start with pseudo count
		for (auto j = 0; j < alphaSize; j++) {
			sum1 += h1[i + j];
			sum2 += h2[i + j];
		}

		for (auto j = 0; j < alphaSize; j++) {
			double m1 = (h1[i + j] + 1.0) / sum1;
			double m2 = (h2[i + j] + 1.0) / sum2;
			double m1AndM2 = m1 + m2;

			d1 += m1 * log(2 * m1 / m1AndM2);
			d2 += m2 * log(2 * m2 / m1AndM2);
		}
	}

	return (d1 + d2) / 2.0;
}

template<class V>
double Statistician<V>::intersectionSimilarity() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		uint64_t s = h1[i] + h2[i];
		// Skip when both entries are zeros.
		if (s != 0) {
			d += 2.0 * std::min(h1[i], h2[i]) / s;
		}
	}
	return d;
}

/**
 * The original statistic is not defined when the two sequences are the same.
 */
template<class V>
double Statistician<V>::kulczynski1Similarity() {
	double d = 0.0;
	double delta = 1.0 / histogramSize;
	for (int i = 0; i < histogramSize; i++) {
		// Skip if both entries are zeros.
		if (h1[i] > 0 || h2[i] > 0) {
			d += (delta + std::min(h1[i], h2[i]))
					/ (delta + absolute(h1[i] - h2[i]));
		}
	}
	return d;
}

/**
 * Note: output of different sequence pair is on the same scale.
 * Comparing any sequence to itself produces the same number.
 */
template<class V>
double Statistician<V>::kulczynski2Similarity() {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d += std::min(h1[i], h2[i]);
	}

	double mu = histogramSize * (mean1 + mean2) / (2 * mean1 * mean2);

	return mu * d;
}

template<class V>
double Statistician<V>::covarianceSimilarityHelper(const V *t1, const V *t2,
		double m1, double m2) {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		d += (t1[i] - m1) * (t2[i] - m2);
	}
	return d / histogramSize;
}

/**
 * Same as covarianceRSimilarity but the mean vector is already calculated.
 */
template<class V>
double Statistician<V>::covarianceSimilarity() {
	return covarianceSimilarityHelper(h1, h2, mean1, mean2);
}

/**
 * Compare the statistic on two histograms to the statistic on
 * the average histogram versus itself.
 */
template<class V>
double Statistician<V>::covarianceRSimilarity() {
	double meanOverall = mean(mean1And2);
	double n = covarianceSimilarityHelper(h1, h2, mean1, mean2);
	double d = covarianceSimilarityHelper(mean1And2, mean1And2, meanOverall,
			meanOverall);

	if (Util::isEqual(d, 0.0)) {
		std::cerr << "Statistician warning at covarianceRSimilarity. ";
		std::cerr << "A sequence is too short. Similarity is assigned zero.";
		std::cerr << std::endl;
		n = 0.0;
	} else {
		n /= d;
	}

	return n;
}

/**
 * Helper to harmonicMeanSimilarity and harmonicMeanRSimilarity
 */
template<class V>
double Statistician<V>::harmonicMeanSimilarityHelper(const V *t1, const V *t2) {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		// Skip row if both entries are zeros.
		if (t1[i] > 0 || t2[i] > 0) {
			d += (t1[i] * t2[i]) / (t1[i] + t2[i]);
		}
	}
	return 2 * d;
}

template<class V>
double Statistician<V>::harmonicMeanSimilarity() {
	return harmonicMeanSimilarityHelper(h1, h2);
}

/**
 * Compare the harmonic mean of two histograms to the harmonic mean of
 * the average histogram versus itself.
 */
template<class V>
double Statistician<V>::harmonicMeanRSimilarity() {
	double n = harmonicMeanSimilarityHelper(h1, h2);
	double d = harmonicMeanSimilarityHelper(mean1And2, mean1And2);

	if (Util::isEqual(d, 0.0)) {
		std::cerr << "Statistician warning at harmonicMeanRSimilarity. ";
		std::cerr << "A sequence is too short. Similarity is assigned zero.";
		std::cerr << std::endl;
		n = 0.0;
	} else {
		n /= d;
	}

	return n;
}

template<class V>
double Statistician<V>::simRatioSimilarity() {
	double dot = 0.0;
	double norm = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		dot += h1[i] * h2[i];
		V diff = h1[i] - h2[i];
		norm += diff * diff;
	}
	double d = dot + sqrt(norm);

	if (Util::isEqual(d, 0.0)) {
		std::cerr << "Error at Sim Ratio. ";
		std::cerr << "The denominator is zero." << std::endl;
		throw std::exception();
	}

	return dot / d;
}

/**
 * Should work on DNA and Proteins
 * t1: histogram of the first sequence
 * t2: histogram of the second sequence
 */
template<class V>
double Statistician<V>::markovSimilarityHelper(const V *t1, const V *t2) {
	double total = 0.0;

	const int alphaSize = Parameters::getAlphabetSize();

	for (auto i = 0; i < histogramSize; i += alphaSize) {
		uint64_t sum1 = alphaSize, sum2 = alphaSize;
		for (auto j = 0; j < alphaSize; j++) {
			sum1 += t1[i + j];
			sum2 += t2[i + j];
		}
		double lsum1 = log(sum1);
		double lsum2 = log(sum2);
		for (auto j = 0; j < alphaSize; j++) {
			// t2 under t1's model
			total += t2[i + j] * (log(t1[i + j] + 1) - lsum1);
			// t1 under t2's model
			total += t1[i + j] * (log(t2[i + j] + 1) - lsum2);
		}
	}

	return total / 2.0;
}

template<class V>
double Statistician<V>::markovSimilarity() {
	return markovSimilarityHelper(h1, h2);
}

/**
 * A version of Markov Similarity that guarantees the similarity between
 * any sequence and itself is the same.
 */
template<class V>
double Statistician<V>::markovRSimilarity() {
	// We subtract (do not divide) here because this statistic is in log scale.
	return markovSimilarityHelper(h1, h2)
			- 0.5
					* (markovSimilarityHelper(h1, h1)
							+ markovSimilarityHelper(h2, h2));
}

/**
 * Should work on DNA and Proteins
 */
template<class V>
double Statistician<V>::simMMSimilarity() {
	double oneUnderOne = 0.0;
	double oneUnderTwo = 0.0;
	double twoUnderOne = 0.0;
	double twoUnderTwo = 0.0;

	const int alphaSize = Parameters::getAlphabetSize();

	for (auto i = 0; i < histogramSize; i += alphaSize) {
		uint64_t sum1 = alphaSize, sum2 = alphaSize;
		for (auto j = 0; j < alphaSize; j++) {
			sum1 += h1[i + j];
			sum2 += h2[i + j];
		}
		double lsum1 = log(sum1);
		double lsum2 = log(sum2);

		for (auto j = 0; j < alphaSize; j++) {
			double hani1 = (log(h1[i + j] + 1) - lsum1);
			double hani2 = (log(h2[i + j] + 1) - lsum2);

			// h1 under h1's model
			oneUnderOne += h1[i + j] * hani1;
			// h1 under h2's model
			oneUnderTwo += h1[i + j] * hani2;
			// h2 under h1's model
			twoUnderOne += h2[i + j] * hani1;
			// h2 under h2's model
			twoUnderTwo += h2[i + j] * hani2;
		}
	}

	uint64_t l1 = sum(h1);
	uint64_t l2 = sum(h2);
	double r = (1.0 / l2) * log(twoUnderOne / twoUnderTwo);
	r += (1.0 / l1) * log(oneUnderTwo / oneUnderOne);
	r /= 2.0;
	return 1.0 - exp(r);
}

/**
 * Calculates the length ratio
 */
template<class V>
double Statistician<V>::lengthRatioSimilarity() {
	uint64_t l1 = sum(h1) + k - 1;
	uint64_t l2 = sum(h2) + k - 1;

	return std::min(l1, l2) / (double) std::max(l1, l2);
}

/**
 * s: a vector of standard deviation of all samples in the data set
 */
template<class V>
double Statistician<V>::seuclideanDistance(const double *s) {
	double d = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		double temp = h1[i] - h2[i];
		if (!Util::isEqual(s[i], 0.0)) {
			d += (temp * temp) / (s[i] * s[i]);
		}
	}
	return sqrt(d);
}

/**
 * Should work on DNA and Protein sequences
 *
 * t1: histogram of the first sequence
 * t2: histogram of the second sequence
 * n: background model, e.g. n[4] = {0.25, 0.25, 0.25, 0.25}
 * keyList: key list in digit format, e.g. getKeysDigitFormat defined in KmerHistogram
 */
template<class V>
double Statistician<V>::d2sSimilarityHelper(const V *t1, const V *t2) {
// Adjust count and calculate statistic
	uint64_t l1 = sum(t1);
	uint64_t l2 = sum(t2);
	if (l1 == 0 || l2 == 0) {
		std::cerr << "Error at d2sSimilarityHelper. ";
		std::cerr << "Sum 1 (l1) or sum 2 (l2) is zero." << std::endl;
		throw std::exception();
	}

	double d2 = 0.0;
	for (int i = 0; i < histogramSize; i++) {
		// Calculate expected values for a word according to the two sequences.
		double e1 = l1;
		double e2 = l2;
		for (int j = 0; j < k; j++) {
			e1 *= background[keyList[i * k + j]];
			e2 *= background[keyList[i * k + j]];
		}
		// Adjust original counts by subtracting the expected values.
		double a1 = t1[i] - e1;
		double a2 = t2[i] - e2;

		double denom = sqrt(a1 * a1 + a2 * a2);
		// Skip undefined rows when both adjusted counts are zeros.
		if (!Util::isEqual(denom, 0.0)) {
			d2 += (a1 * a2) / denom;
		} else {
			std::cout << "Skipped a row" << std::endl;
		}
	}

	return d2;
}

template<class V>
double Statistician<V>::d2sSimilarity() {
	return d2sSimilarityHelper(h1, h2);
}

/**
 * Compare the statistic on two histograms to the statistic on
 * the average histogram versus itself.
 */
template<class V>
double Statistician<V>::d2sRSimilarity() {
	return d2sSimilarityHelper(h1, h2)
			/ d2sSimilarityHelper(mean1And2, mean1And2);
}

/**
 * Should work on DNA and Protein sequences
 */
template<class V>
double Statistician<V>::d2starSimilarity() {
	// Calculate probability vectors based on monomers.
	const int alphaSize = Parameters::getAlphabetSize();
	uint64_t s1 = sum(mono1, alphaSize);
	uint64_t s2 = sum(mono2, alphaSize);

	if (s1 == 0 || s2 == 0) {
		std::cerr << "Error at d2starSimilarity. ";
		std::cerr << "Sum 1 (s1) or sum 2 (s2) is zero." << std::endl;
		throw std::exception();
	}

	// Calculate probability vector from two sequences.
	double p[alphaSize];
	uint64_t s = s1 + s2;
	for (int i = 0; i < alphaSize; i++) {
		p[i] = ((double) mono1[i] + mono2[i] + 1.0) / (s + alphaSize);
	}

	// Adjust count and calculate statistic
	uint64_t l1 = sum(h1);
	uint64_t l2 = sum(h2);
	if (l1 == 0 || l2 == 0) {
		std::cerr << "Error at d2sSimilarity. ";
		std::cerr << "Sum 1 (l1) or sum 2 (l2) is zero." << std::endl;
		throw std::exception();
	}

	double d2 = 0.0;
	double l = sqrt(l1 * l2);
	for (int i = 0; i < histogramSize; i++) {
		double e1 = l1;
		double e2 = l2;
		double e = l;
		for (int j = 0; j < k; j++) {
			// Calculate expected values for a word according to the background model.
			e1 *= background[keyList[i * k + j]];
			e2 *= background[keyList[i * k + j]];
			// Calculate expected values for a word according to the two sequences.
			e *= p[keyList[i * k + j]];
		}
		// Adjust original counts by subtracting the expected values.
		double a1 = h1[i] - e1;
		double a2 = h2[i] - e2;

		// Skip undefined rows when both adjusted counts are zeros.
		if (!Util::isEqual(e, 0.0)) {
			d2 += ((a1 * a2) / e);
		} else {
			std::cout << "Skipped a row" << std::endl;
		}
	}
	return d2;
}

/**
 * r: a vector of the resulting statistics
 */
template<class V>
void Statistician<V>::calculateAll(std::vector<double> &r) {
	for (int i = 0; i < Stat::Dist_NUM; i++) {
		r.push_back((this->*methodList[i])());
	}

	for (int i = Stat::Dist_NUM + 1; i < Stat::ALL_NUM; i++) {
		r.push_back((this->*methodList[i])());
	}
}

/**
 * s: a vector of method indexes
 * r: a vector of the resulting statistics
 */
template<class V>
void Statistician<V>::calculate(std::vector<int> &s, std::vector<double> &r) {
	for (int i : s) {
		if ((i >= 0 && i < Stat::Dist_NUM)
				|| (i > Stat::Dist_NUM && i < Stat::ALL_NUM)) {
			r.push_back((this->*methodList[i])());
		} else {
			std::cerr << "Statistician error: Invalid statistic index.";
			std::cerr << std::endl;
			throw std::exception();
		}
	}
}

template<class V>
void Statistician<V>::calculate(const int *s, int size, double *r) {
	for (int i = 0; i < size; i++) {
		r[i] = (this->*methodList[s[i]])();
	}
}
