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
 * matrix.h
 *
 * Created on: May 10, 2017
 * Modified on: 4/14/2020
 * Author: Hani Z. Girgis, The Bioinformatics Toolsmith Laboratory
 * Author: Robert Geraghty, The Bioinformatics Toolsmith Laboratory
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <string>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <climits>

#include "Util.h"

class Matrix {
private:
	double *m;
	int numRow;
	int numCol;
	inline int indexOf(int i, int j) const {
		return i * numCol + j;
	}

public:
	Matrix();
	Matrix(int r, int c, double init = 0);
	Matrix(const Matrix&);
	Matrix(Matrix&&);

	~Matrix();

	Matrix& operator=(const Matrix &n);

	Matrix operator+(const Matrix &n) const;
	Matrix operator-(const Matrix &n) const;
	Matrix operator*(const Matrix &n) const;
	Matrix operator~() const; // Transpose
	Matrix operator!() const; // Inverse
	Matrix pseudoInverse() const;

	// Returns a matrix allocated on the heap
	Matrix* times(const Matrix*) const;
	Matrix* transpose() const;
	Matrix* inverse() const;
	Matrix* pseudoInverseHeap() const;

	// Make a sub-matrix efficiently
	Matrix* subMatrix(const int[], int) const;
	Matrix* subMatrixByColHeap(const int[], int) const;
	Matrix subMatrixByRow(const int[], int) const;
	Matrix subMatrixByCol(const int[], int) const;

	// Helper to pseudoInverse() and pseudoInverseHeap()
	Matrix inverseHelper() const;

	/**
	 * Setters and getters
	 */
	// Credit: https://www.learncpp.com/cpp-tutorial/99-overloading-the-parenthesis-operator/
	double& operator()(int row, int col);
	// Credit: https://www.learncpp.com/cpp-tutorial/99-overloading-the-parenthesis-operator/
	// for constant objects
	const double& operator()(int row, int col) const;
	// Same function as the () operator
	double& item(int, int);
	// Same function as the () operator
	const double& item(int row, int col) const;
	// Access item after boundary check
	double& at(int, int);
	// Access item after boundary check
	const double& at(int row, int col) const;

	double* getArray();
	const double* getArray() const;

	std::vector<double> getRow(int r) const;
	void setRow(int, std::vector<double>&);
	void setRow(int, int, std::vector<double>&);
	void setCol(int, std::vector<double>&);
	int getNumRow() const;
	int getNumCol() const;

	// Print methods
	void print();
	void printToFile(std::string);

	// Fill methods
	void userFill();
	void randFill(double low, double high);
	void fileFill(std::string filename);

	std::vector<int> getNonZerosRows();

	/**
	 * Append ones column
	 */
	Matrix appendOnesColumn();

	void clearData();
};

#endif /* MATRIX_H_ */
