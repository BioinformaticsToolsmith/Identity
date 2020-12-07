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
 * BestFirst.cpp
 *
 *  Created on: Apr 28, 2020
 *      Author: Dr. Hani Z. Girgis
 *
 */

/**
 * isHigher: true for classification (accuracy) and false for regression (mse or mae)
 */
template<class T>
BestFirst<T>::BestFirst(const Matrix &f, const Matrix &l,
		std::vector<Feature*> &h, MakeTransformer f1, Evaluate f2,
		IsNewBetter f3, bool isHigher, int threadNum, int minFeatIn,
		double lowest, int k) {
	fList = copy(h);
	fNum = f.getNumCol();
	makeTransformer = f1;
	evaluate = f2;
	isNewBetter = f3;
	limit = lowest;
	minFeat = minFeatIn;

	isHigherBetter = isHigher;

	Node empty;
	best = empty;
	maximum = lowest; // This is ok for the current evaluation measures
	open[empty] = lowest;

	/**
	 * The search algorithm
	 * See page 293 of R. Kohavi, G.H. John/Artificial Intelligence 97 (1997) 273-324
	 *
	 */
	int eCount = 0;
	while ((best.getSize() < minFeat || eCount < k) && !open.empty()) {
		// Step 2
		std::pair<Node, double> v = findOptimum();

		// Step 3
		open.erase(v.first);
		closed[v.first] = v.second;
		// Step 4
		if (isNewBetter(v.second, maximum)
				|| (best.getSize() < minFeat && eCount >= k)) {
			best = v.first;
			maximum = v.second;
			eCount = 0;

			// Print new features
			std::cout << "Better performance of: " << maximum << std::endl;
			const int *list = best.getList();
			int s = best.getSize();
			for (int a = 0; a < s; a++) {
				std::cout << "\t" << fList.at(list[a])->getName() << std::endl;
			}
		}
		// Step 5
		std::vector < Node > eList = v.first.expand(fNum);
		eCount++;

		// Step 6
		int childNum = eList.size();
		// double resultList[childNum] { 0.0 };
		std::vector<double> resultList(childNum, 0.0);
#pragma omp parallel for schedule(static) num_threads(threadNum)
		for (int i = 0; i < childNum; i++) {
			Node child = eList[i];

			if (open.find(child) == open.end()
					&& closed.find(child) == closed.end()) {
				Matrix t =
						f.subMatrixByCol(child.getList(), child.getSize()).appendOnesColumn();
				try {
					T predictor = makeTransformer(t, l);
					double e = evaluate(l, predictor.transform(t));
					resultList[i] = e;
				} catch (const std::exception &e) {
#pragma omp critical
					{
						std::cerr << "A feature set was ignored:" << std::endl;
						const int *list = child.getList();
						for (int a = 0; a < child.getSize(); a++) {
							std::cout << "\t" << list[a] << ": ";
							std::cout << fList.at(list[a])->getName()
									<< std::endl;
						}

						std::cerr
								<< "The problematic matrix is written to exception.txt";
						std::cerr << std::endl;
						t.printToFile("exception.txt");

						throw std::exception();
					}
				}
			}
		}

		for (int i = 0; i < childNum; i++) {
			Node child = eList[i];
			if (open.find(child) == open.end()
					&& closed.find(child) == closed.end()) {
				open[child] = resultList[i];
			}
		}

	}

	std::cout << "Selected statistics:" << std::endl;
	// Print new features
	const int *list = best.getList();
	int s = best.getSize();
	for (int a = 0; a < s; a++) {
		std::cout << "\t" << fList.at(list[a])->getName() << std::endl;
	}

	// Empty open and closed
	open.clear();
	closed.clear();
}

/**
 * Constructor of optimized features
 */
template<class T>
BestFirst<T>::BestFirst(std::vector<Feature*> &h) {
	fList = copy(h);

	// Count selected features
	int s = 0;
	for (auto f : fList) {
		if (f->getIsSelected()) {
			s++;
		}
	}

	// Construct best node
	int *l = new int[s];
	int i = 0;
	for (auto f : fList) {
		if (f->getIsSelected()) {
			l[i] = f->getTableIndex();
			i++;
		}
	}

	best = Node(l, s);
}

template<class T>
BestFirst<T>::~BestFirst() {
	for (auto ptr : fList) {
		delete ptr;
	}
	fList.clear();
}

template<class T>
std::pair<Node, double> BestFirst<T>::findOptimum() {
	if (open.empty()) {
		std::cerr << "BestFirst error: " << std::endl;
		std::cerr << "Cannot find optimum on empty Open set." << std::endl;
		throw std::exception();
	}

	auto it = open.begin();
	Node k;
	double v = limit;
	while (it != open.end()) {
		if ((isHigherBetter && it->second > v)
				|| (!isHigherBetter && it->second < v)) {
			k = it->first;
			v = it->second;
		}
		it++;
	}
	return std::make_pair(k, v);
}

template<class T>
Matrix BestFirst<T>::transform(const Matrix &m) {
	const int *l = best.getList();
	int s = best.getSize();
	return m.subMatrixByCol(l, s).appendOnesColumn();
}

template<class T>
std::vector<Feature*> BestFirst<T>::getFeatureList() {
	const int *l = best.getList();
	int s = best.getSize();
	std::vector<Feature*> temp = copy(fList);
	// Mark selected (and automatically needed) features
	for (int i = 0; i < s; i++) {
		temp.at(l[i])->setIsSelected();
	}

	std::vector<Feature*> r;
	for (auto ptr : temp) {
		if (ptr->getIsSelected() || ptr->getIsNeeded()) {
			r.push_back(ptr);
		} else {
			delete ptr;
		}
	}
	r.shrink_to_fit();

	return r;
}

template<class T>
const Node& BestFirst<T>::getBest() const {
	return best;
}
