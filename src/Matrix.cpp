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
 * matrix.cpp
 *
 * Created on: May 10, 2017
 * Rewritten on: April 15, 2020 By Dr. Hani Z. Girgis, The Bioinformatics Toolsmith Laboratory
 * Author: Robert Geraghty, The Bioinformatics Toolsmith Laboratory
 *
 * ToDo:
 * [OK] Test the inverse on this matrix https://www.mathsisfun.com/algebra/matrix-inverse-row-operations-gauss-jordan.html
 * [OK] Review gaussJordanInverse()
 * [OK] Add empty constructor
 *
 * [OK] Make sure that there is no memory leak, specially when we construct a matrix and return it from a method.
 * [OK] Stop using indexOf because it depends on the column size of a matrix.
 * [OK] Define another method similar to at() but does not check boundaries.
 *
 *
 * [  ] Make it a template class using float? Not sure if it wall affect precision
 */

#include "Matrix.h"

Matrix::Matrix() :
		numRow(0), numCol(0) {
	m = nullptr;
}

Matrix::Matrix(int r, int c, double init) :
		numRow(r), numCol(c) {
	int size = r * c;
	m = new double[size];
	for (int i = 0; i < size; i++) {
		m[i] = init;
	}
}

Matrix::Matrix(const Matrix &o) :
		numRow(o.getNumRow()), numCol(o.getNumCol()) {
	int size = numRow * numCol;
	m = new double[size];
	const double *n = o.getArray();
	for (int i = 0; i < size; i++) {
		m[i] = n[i];
	}
}

Matrix& Matrix::operator=(const Matrix &o) {
	numRow = o.getNumRow();
	numCol = o.getNumCol();
	int size = numRow * numCol;
	m = new double[size];
	const double *n = o.getArray();
	for (int i = 0; i < size; i++) {
		m[i] = n[i];
	}

	return *this;
}

Matrix::Matrix(Matrix &&o) :
		numRow(o.getNumRow()), numCol(o.getNumCol()) {
	m = o.m;
	o.m = nullptr;
}

Matrix::~Matrix() {
	if (m != nullptr) {
		delete[] m;
		m = nullptr;
	}
}

double& Matrix::operator()(int r, int c) {
	return m[indexOf(r, c)];
}

double& Matrix::item(int r, int c) {
	return m[indexOf(r, c)];
}

// Same as the () operator, but it checks the bounds.
double& Matrix::at(int r, int c) {
	if (r < 0 || r >= numRow) {
		std::cerr << "Matrix error: Row index is out of range." << std::endl;
		throw std::exception();
	}
	if (c < 0 || c >= numCol) {
		std::cerr << "Matrix error: Column index is out of range." << std::endl;
		std::cerr << "Number of column is " << numCol << ", but received " << c;
		std::cerr << std::endl;
		throw std::exception();
	}

	return m[indexOf(r, c)];
}

const double& Matrix::operator()(int r, int c) const {
	return m[indexOf(r, c)];
}

const double& Matrix::item(int r, int c) const {
	return m[indexOf(r, c)];
}

// Same as the () operator, but it checks the bounds.
const double& Matrix::at(int r, int c) const {
	if (r < 0 || r >= numRow) {
		std::cerr << "Matrix error: Row index is out of range." << std::endl;
		throw std::exception();
	}
	if (c < 0 || c >= numCol) {
		std::cerr << "Matrix error: Column index is out of range." << std::endl;
		throw std::exception();
	}

	return m[indexOf(r, c)];
}

std::vector<double> Matrix::getRow(int r) const {
	double *rowStart = &m[indexOf(r, 0)];
	std::vector<double> row(rowStart, rowStart + numCol);
	return row;
}

/**
 * This method assign a whole row
 */
void Matrix::setRow(int r, std::vector<double> &valList) {
	if (valList.size() != numCol) {
		std::cerr << "The size of the row does not match that of the matrix.";
		std::cerr << std::endl;
		throw std::exception();
	}
	int r0 = indexOf(r, 0);
	for (int c = 0; c < numCol; c++) {
		m[r0 + c] = valList[c];
	}
}

/**
 * Make a sub-matrix by row efficiently
 * indexList: array contained desired indexes
 * indexSize: size of indexList
 * Memory: It is a client's responsibility to free memory allocated to the sub-matrix
 */
Matrix* Matrix::subMatrix(const int indexList[], int indexSize) const {
	Matrix *result = new Matrix(indexSize, numCol);
	double *n = result->getArray();
	for (int i = 0; i < indexSize; i++) {
		int r = indexList[i];
		if (r < numRow) {
			for (int c = 0; c < numCol; c++) {
				n[indexOf(i, c)] = m[indexOf(r, c)];
			}
		} else {
			std::cerr
					<< "Matrix error: Row index in sub-matrix is out of range.";
			std::cerr << std::endl;
			throw std::exception();
		}
	}

	return result;
}

/**
 * Make a sub-matrix by row efficiently
 * indexList: array contained desired indexes
 * indexSize: size of indexList
 */
Matrix Matrix::subMatrixByRow(const int indexList[], int indexSize) const {
	Matrix result(indexSize, numCol);
	double *n = result.getArray();
	for (int i = 0; i < indexSize; i++) {
		int r = indexList[i];
		if (r < numRow) {
			// ToDo: Calculate the first ic and rc and then add one to each
			for (int c = 0; c < numCol; c++) {
				n[indexOf(i, c)] = m[indexOf(r, c)];
			}
		} else {
			std::cerr
					<< "Matrix error: Row index in sub-matrix is out of range.";
			std::cerr << std::endl;
			throw std::exception();
		}
	}

	return result;
}

/**
 * Make a sub-matrix by column efficiently
 * indexList: array contained desired column indexes
 * indexSize: size of indexList
 */
Matrix Matrix::subMatrixByCol(const int indexList[], int indexSize) const {
	Matrix result(numRow, indexSize);
	double *n = result.getArray();
	for (int i = 0; i < indexSize; i++) {
		int c = indexList[i];
		if (c < numCol) {
			int ri = i;
			int rc = c;
			for (int r = 0; r < numRow; r++) {
				n[ri] = m[rc];
				ri += indexSize; // Next entry in a column of new matrix
				rc += numCol; // Next entry in a column of original matrix
			}
		} else {
			std::cerr << "Matrix error: ";
			std::cerr << "Column index in sub-matrix is out of range.";
			std::cerr << std::endl;
			throw std::exception();
		}
	}

	return result;
}

Matrix* Matrix::subMatrixByColHeap(const int indexList[], int indexSize) const {
	Matrix *result = new Matrix(numRow, indexSize);
	double *n = result->getArray();
	for (int i = 0; i < indexSize; i++) {
		int c = indexList[i];
		if (c < numCol) {
			int ri = i;
			int rc = c;
			for (int r = 0; r < numRow; r++) {
				n[ri] = m[rc];
				ri += indexSize; // Next entry in a column of new matrix
				rc += numCol; // Next entry in a column of original matrix
			}
		} else {
			std::cerr << "Matrix error: ";
			std::cerr << "Column index in sub-matrix is out of range.";
			std::cerr << std::endl;
			throw std::exception();
		}
	}

	return result;
}

/**
 * This method assign a part of a row
 * r: row number
 * f: starting column number
 */
void Matrix::setRow(int r, int f, std::vector<double> &valList) {
	if (f + valList.size() > numCol) {
		std::cerr << "Matrix Error: Exceeding row size of the matrix.";
		std::cerr << std::endl;
		throw std::exception();
	}
	int rf = indexOf(r, f);
	for (int c = 0; c < valList.size(); c++) {
		m[rf + c] = valList[c];
	}
}

/**
 * This method assign a whole column
 */
void Matrix::setCol(int c, std::vector<double> &valList) {
	if (valList.size() != numRow) {
		std::cerr
				<< "Matrix Error: The size of the column does not match that of the matrix.";
		std::cerr << std::endl;
		throw std::exception();
	}

	for (int r = 0; r < numRow; r++) {
		m[indexOf(r, c)] = valList[r];
	}
}

int Matrix::getNumRow() const {
	return numRow;
}

int Matrix::getNumCol() const {
	return numCol;
}

double* Matrix::getArray() {
	if (m == nullptr) {
		std::cerr << "Matrix error: Array is null." << std::endl;
		throw std::exception();
	}
	return m;
}

const double* Matrix::getArray() const {
	if (m == nullptr) {
		std::cerr << "Matrix error: Array is null." << std::endl;
		throw std::exception();
	}
	return m;
}

Matrix Matrix::operator+(const Matrix &other) const {
	if (numCol != other.numCol || numRow != other.numRow) {
		std::cerr << "Invalid input: Matrix dimension mismatch." << std::endl;
		throw std::exception();
	}

	Matrix result(numRow, numCol);
	const double *n = other.getArray();
	double *r = result.getArray();
	int total = numRow * numCol;
	for (int i = 0; i < total; i++) {
		r[i] = m[i] + n[i];
	}
	return result;
}

Matrix Matrix::operator-(const Matrix &other) const {
	if (numCol != other.numCol || numRow != other.numRow) {
		std::cerr << "Invalid input: Matrix dimension mismatch." << std::endl;
		throw std::exception();
	}

	Matrix result(numRow, numCol);
	const double *n = other.getArray();
	double *r = result.getArray();
	int total = numRow * numCol;
	for (int i = 0; i < total; i++) {
		r[i] = m[i] - n[i];
	}
	return result;
}

/**
 * Returns an object allocated on the heap
 * It is the client's responsibility to free memory allocated to the result
 */
Matrix* Matrix::times(const Matrix *other) const {
	if (numCol != other->getNumRow()) {
		std::cerr << "Matrix error: ";
		std::cerr << "In times because matrix dimension mismatch." << std::endl;
		std::cerr << "Columns of the first matrix is: " << numCol << std::endl;
		std::cerr << "Columns of the second matrix is: " << other->getNumRow();
		std::cerr << std::endl;
		throw std::exception();
	}

	Matrix *result = new Matrix(numRow, other->getNumCol());

	double curSum = 0.0;
	for (int i = 0; i < result->getNumRow(); i++) {
		for (int j = 0; j < result->getNumCol(); j++) {
			curSum = 0;
			for (int k = 0; k < numCol; k++) {
				curSum = curSum + item(i, k) * other->item(k, j);
			}
			result->item(i, j) = curSum;
		}
	}

	return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
	if (numCol != other.getNumRow()) {
		std::cerr << "Invalid input: Matrix dimension mismatch." << std::endl;
		throw std::exception();
	}

	double curSum = 0.0;
	Matrix result(numRow, other.getNumCol());

	for (int i = 0; i < result.getNumRow(); i++) {
		for (int j = 0; j < result.getNumCol(); j++) {
			curSum = 0;
			for (int k = 0; k < numCol; k++) {
				curSum = curSum + m[indexOf(i, k)] * other(k, j);
			}
			result(i, j) = curSum;
		}
	}

	return result;
}

/**
 * Returns an object allocated on the heap
 * It is the client's responsibility to free memory allocated to it
 */
Matrix* Matrix::transpose() const {
	Matrix *t = new Matrix(numCol, numRow);

	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			t->item(j, i) = item(i, j);
		}
	}

	return t;
}

/**
 * Transpose
 */
Matrix Matrix::operator~() const {
	Matrix t = Matrix(numCol, numRow);
	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			t(j, i) = item(i, j);
		}
	}
	return t;
}

/**
 * Algorithm: https://online.stat.psu.edu/statprogram/reviews/matrix-algebra/gauss-jordan-elimination
 */
Matrix Matrix::inverseHelper() const {
	if (numRow != numCol) {
		std::cerr << "Matrix error: ";
		std::cerr << "Cannot take the inverse of a non-squared matrix.";
		std::cerr << std::endl;
		throw std::exception();
	}

	// m -> m ID
	Matrix matrix = Matrix(numRow, 2 * numCol);
	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			matrix(i, j) = item(i, j);
			if (i == j) {
				matrix(i, j + numCol) = 1.0;
			}
		}
	}

	// Add/subtract multiples of the top row to the other rows so that
	// all other entries in the column containing the top row's leading entry are all zero.
	// At this point, pivots (leading entries) are not zeros.
	for (int i = 0; i < numRow; i++) {
		// Handle the case when the pivot is zero
		// Find another row with non-zero in the same column as the pivot
		// and add this row to the current row
		if (Util::isEqual(matrix(i, i), 0.0)) {
			for (int j = 0; j < numRow; j++) {
				if (i != j && !Util::isEqual(matrix(j, i), 0.0)) {
					for (int k = 0; k < 2 * numRow; k++) {
						matrix(i, k) += matrix(j, k);
					}
					break;
				}
			}
		}

		for (int j = 0; j < numRow; j++) {
			if (j != i) {
				if (!Util::isEqual(matrix(j, i), 0.0)) {
					double temp = matrix(j, i) / matrix(i, i);

					for (int k = 0; k < 2 * numRow; k++) {
						if (k == i) {
							matrix(j, k) = 0.0;
						} else {
							matrix(j, k) -= (matrix(i, k) * temp);
						}
					}
				}
			}
		}
	}

	// Divide row element by the diagonal element (pivot).
	// At this point, pivots are not zeros.
	for (int i = 0; i < numRow; i++) {
		double temp = matrix(i, i);

		if (!Util::isEqual(temp, 1.0)) {
			for (int j = 0; j < 2 * numRow; j++) {
				matrix(i, j) /= temp;
			}
		}
	}

	// Now check to make sure the original matrix is an identity matrix.
	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			if (i == j && !Util::isEqual(matrix(i, j), 1.0)) {
				std::cerr << "Matrix error: No inverse found." << std::endl;
				std::cerr << "\tExpected 1.0 but found: " << matrix(i, j)
						<< std::endl;
				throw std::exception();
			}

			if (i != j && !Util::isEqual(matrix(i, j), 0.0)) {
				std::cerr << "Matrix error: No inverse found." << std::endl;
				std::cerr << "\tExpected 0.0 but found: " << matrix(i, j)
						<< std::endl;
				throw std::exception();
			}
		}
	}

	return matrix;
}

/**
 * Gauss-Jordan inverse
 */
Matrix Matrix::operator!() const {
	if (numRow != numCol) {
		std::cerr << "Matrix error: Not a squared matrix." << std::endl;
		throw std::exception();
	}

	int indexList[numCol];
	for (int i = 0; i < numCol; i++) {
		indexList[i] = numCol + i;
	}

	return inverseHelper().subMatrixByCol(indexList, numCol);
}

/**
 * The returned matrix is allocated on the heap
 * It is the client's responsibility to free memory allocated to the result
 */
Matrix* Matrix::inverse() const {
	if (numRow != numCol) {
		std::cerr << "Matrix error: Not a squared matrix." << std::endl;
		throw std::exception();
	}

	int indexList[numCol];
	for (int i = 0; i < numCol; i++) {
		indexList[i] = numCol + i;
	}

	return inverseHelper().subMatrixByColHeap(indexList, numCol);

}

Matrix Matrix::pseudoInverse() const {
	if (numRow >= numCol) {
		Matrix trans = this->operator~();
		return !(trans * *this) * trans;
	} else {
		Matrix trans = this->operator~();
		Matrix psuedoInv = trans * (!this->operator*(trans));
		return psuedoInv;
	}
}

Matrix* Matrix::pseudoInverseHeap() const {
	if (numRow >= numCol) {
		Matrix *temp = transpose();
		Matrix *transByOrig = temp->times(this);
		Matrix *transByOrigInverse = transByOrig->inverse();
		Matrix *psuedoInv = transByOrigInverse->times(temp);

		delete temp;
		delete transByOrig;
		delete transByOrigInverse;

		return psuedoInv;
	} else {
		Matrix *trans = transpose();
		Matrix *origByTrans = this->times(trans);
		Matrix *origByTransInverse = origByTrans->inverse();
		Matrix *psuedoInv = trans->times(origByTransInverse);

		delete trans;
		delete origByTrans;
		delete origByTransInverse;

		return psuedoInv;
	}
}

void Matrix::print() {
	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			std::cout << m[indexOf(i, j)] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

/**
 * Print in Matlab format
 */
void Matrix::printToFile(std::string fileName) {
	std::ofstream out(fileName.c_str());

	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			out << std::right << std::fixed;
			out << std::setprecision(16) << std::setw(7) << m[indexOf(i, j)];
			if (j < numCol - 1) {
				out << "\t";
			} else {
				out << std::endl;
			}
		}
	}

	out.close();
}

void Matrix::randFill(double low, double high) {
	double x;
	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			x = ((double) rand() * (high - low)) / (double) RAND_MAX + low;
			m[indexOf(i, j)] = x;
		}
	}
}

void Matrix::userFill() {
	double val;
	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			std::cout << "input value for cell (" << i << ", " << j << ")?\n";
			std::cin >> val;
			std::cout << std::endl;
			m[indexOf(i, j)] = val;
		}
	}
}

void Matrix::fileFill(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (!infile) {
		std::cerr << "Matrix Error: File read fail" << std::endl;
		throw std::exception();
	}

	std::string line;
	int i = -1;
	while (getline(infile, line)) {
		i++;
		if (i >= numRow) {
			std::cerr << "Matrix Error: ";
			std::cerr << "Matrix in file has incorrect number of rows";
			std::cerr << std::endl;
			throw std::exception();
		}

		double num;
		std::istringstream iss(line);
		int j = -1;
		while (iss >> num) {
			j++;
			if (j >= numCol) {
				std::cerr << "Matrix Error: ";
				std::cerr << "Matrix in file has incorrect number of columns";
				std::cerr << std::endl;
				throw std::exception();
			}
			m[indexOf(i, j)] = num;
		}
		j = 0;
	}
	i = 0;
}

/**
 * Find rows that are not all-zeros
 */
std::vector<int> Matrix::getNonZerosRows() {
	std::vector<int> r;
	r.reserve(numRow);
	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numCol; j++) {
			if (!Util::isEqual(m[indexOf(i, j)], 0.0)) {
				r.push_back(i);
				break;
			}
		}
	}
	return r;
}

/**
 * The ones column is the first column
 */
Matrix Matrix::appendOnesColumn() {
	Matrix h(numRow, numCol + 1, 1.0);
	for (int r = 0; r < numRow; r++) {
		for (int c = 0; c < numCol; c++) {
			h(r, c + 1) = m[indexOf(r, c)];
		}
	}
	return h;
}
