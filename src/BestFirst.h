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
 * BestFirst.h
 *
 *  Created on: Apr 28, 2020
 *      Author: Dr. Hani Zakaria Girgis
 */

#ifndef BESTFIRST_H_
#define BESTFIRST_H_

#include <unordered_map>
#include <vector>

#include "Node.h"
#include "Matrix.h"
#include "ITransformer.h"
#include "Feature.h"
#include "Parameters.h"

template<class T>
class BestFirst: public ITransformer {

	// A function for make a classifier or a regression model
	typedef T (*MakeTransformer)(const Matrix&, const Matrix&);
	// A function to measure performance, e.g. accuracy or squared mean error
	typedef double (*Evaluate)(const Matrix&, const Matrix&);
	// A function to determine if performance is improving
	// !! The order of these two parameters is very important
	typedef bool (*IsNewBetter)(double newV, double oldV);

private:
	std::unordered_map<Node, double, NodeHasher> open;
	std::unordered_map<Node, double, NodeHasher> closed;

	Node best;

	double maximum;
	int fNum;
	std::vector<Feature*> fList;
	double limit;

	MakeTransformer makeTransformer;
	Evaluate evaluate;
	IsNewBetter isNewBetter;

	bool isHigherBetter = true;
	int minFeat = 0;

	/**
	 * Return the best node in the open set depending on the
	 * evaluation measure used and whether it is a maximization
	 * or a minimization process
	 * It return the value associated by this node
	 */
	std::pair<Node, double> findOptimum();

public:
	BestFirst(const Matrix&, const Matrix&, std::vector<Feature*>&,
			MakeTransformer, Evaluate, IsNewBetter, bool, int, int, double lowest = 0.0,
			int stop = 3);
	BestFirst(std::vector<Feature*> &h);
	virtual ~BestFirst();
	virtual Matrix transform(const Matrix&);
	virtual std::vector<Feature*> getFeatureList();
	const Node& getBest() const;
};

#include "BestFirst.cpp"

#endif /* BESTFIRST_H_ */
