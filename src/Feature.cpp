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
 * Feature.cpp
 *
 *  Created on: May 19, 2020
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "Feature.h"

/**
 * n: Name
 * i: Distance (true) or similarity (false)
 * f: Function index used by statistician
 */
Feature::Feature(int f, std::string n, bool i) :
		funIndex(f), name(n), isDistance(i) {
	isNormalized = false;
	normP1 = 0.0;
	normP2 = 0.0;
	tableIndex = -1;
	isSelected = false;
	isNeeded = false;
	isConverted = false;
	w = 0.0;
}

Feature::~Feature() {
}

Feature::Feature(const Feature &other) :
		funIndex(other.funIndex), name(other.name), isDistance(other.isDistance) {
	isNormalized = other.isNormalized;
	normP1 = other.normP1;
	normP2 = other.normP2;
	tableIndex = other.tableIndex;
	isSelected = other.isSelected;
	isNeeded = other.isNeeded;
	isConverted = other.isConverted;
	w = other.w;
}

Feature::Feature(Feature &&other) :
		funIndex(other.funIndex), name(other.name), isDistance(other.isDistance) {
	isNormalized = other.isNormalized;
	normP1 = other.normP1;
	normP2 = other.normP2;
	tableIndex = other.tableIndex;
	isSelected = other.isSelected;
	isNeeded = other.isNeeded;
	isConverted = other.isConverted;
	w = other.w;
}

bool Feature::getIsNormalized() const {
	return isNormalized;
}

void Feature::setIsNormalized(bool isNormalized) {
	this->isNormalized = isNormalized;
}

double Feature::getNormP1() const {
	return normP1;
}

void Feature::setNormP1(double normP1) {
	this->normP1 = normP1;
}

double Feature::getNormP2() const {
	return normP2;
}

int Feature::getTableIndex() const {
	return tableIndex;
}

const std::string& Feature::getName() const {
	return name;
}

const int Feature::getFunIndex() const {
	return funIndex;
}

void Feature::setTableIndex(int tableIndex) {
	this->tableIndex = tableIndex;
}

void Feature::setNormP2(double normP2) {
	this->normP2 = normP2;
}

int Feature::getNumOfComp() {
	return 0;
}

bool Feature::getIsDistance() const {
	return isDistance;
}

bool Feature::getIsNeeded() const {
	return isNeeded;
}

void Feature::setIsNeeded() {
	isNeeded = true;
}

bool Feature::getIsSelected() const {
	return isSelected;
}

void Feature::setIsSelected() {
	isSelected = true;
}

int Feature::getCompOneIndex() {
	std::cerr << "Feature error: " << std::endl;
	std::cerr << "A feature does not have components.";
	std::cerr << std::endl;
	throw std::exception();
}

int Feature::getCompTwoIndex() {
	std::cerr << "Feature error: " << std::endl;
	std::cerr << "A feature does not have components.";
	std::cerr << std::endl;
	throw std::exception();
}

void Feature::setCompOne(Feature *f) {
	std::cerr << "Feature error: " << std::endl;
	std::cerr << "A feature does not have components.";
	std::cerr << std::endl;
	throw std::exception();
}

void Feature::setCompTwo(Feature *f) {
	std::cerr << "Feature error: " << std::endl;
	std::cerr << "A feature does not have components.";
	std::cerr << std::endl;
	throw std::exception();
}

bool Feature::getIsConverted() const {
	return isConverted;
}

void Feature::setIsConverted(bool isConverted) {
	this->isConverted = isConverted;
}

double Feature::getW() const {
	return w;
}

void Feature::setW(double w) {
	this->w = w;
}
