/*
 MeShClust v3.0 clusters sequences using the mean shift algorithm and alignment-free identity scores.

 Copyright (C) 2020-2022 Hani Z. Girgis, PhD

 Academic use: Affero General Public License version 1.

 Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 */

/*
 * MeShClust.cpp
 *
 *  Created on: December 23, 2020
 *      Author: Hani Z. Girgis, PhD
 *
 */

#include <iostream>
#include <thread>
#include <algorithm>

#include "GMM.h"
#include "MeanShiftLarge.h"
#include "../Util.h"
#include "../FastaReader.h"
#include "../IdentityCalculator.h"
#include "../IdentityCalculator1.h"
#include "../Parameters.h"

const char *agplv1 =
		R"(AFFERO GENERAL PUBLIC LICENSE
Version 1, March 2002

Copyright © 2002 Affero Inc.
510 Third Street - Suite 225, San Francisco, CA 94107, USA

This license is a modified version of the GNU General Public License copyright (C) 1989, 1991 Free Software Foundation, Inc. made with their permission. Section 2(d) has been added to cover use of software over a computer network.

Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed.

Preamble

The licenses for most software are designed to take away your freedom to share and change it. By contrast, the Affero General Public License is intended to guarantee your freedom to share and change free software--to make sure the software is free for all its users. This Public License applies to most of Affero's software and to any other program whose authors commit to using it. (Some other Affero software is covered by the GNU Library General Public License instead.) You can apply it to your programs, too.

When we speak of free software, we are referring to freedom, not price. This General Public License is designed to make sure that you have the freedom to distribute copies of free software (and charge for this service if you wish), that you receive source code or can get it if you want it, that you can change the software or use pieces of it in new free programs; and that you know you can do these things.

To protect your rights, we need to make restrictions that forbid anyone to deny you these rights or to ask you to surrender the rights. These restrictions translate to certain responsibilities for you if you distribute copies of the software, or if you modify it.

For example, if you distribute copies of such a program, whether gratis or for a fee, you must give the recipients all the rights that you have. You must make sure that they, too, receive or can get the source code. And you must show them these terms so they know their rights.

We protect your rights with two steps: (1) copyright the software, and (2) offer you this license which gives you legal permission to copy, distribute and/or modify the software.

Also, for each author's protection and ours, we want to make certain that everyone understands that there is no warranty for this free software. If the software is modified by someone else and passed on, we want its recipients to know that what they have is not the original, so that any problems introduced by others will not reflect on the original authors' reputations.

Finally, any free program is threatened constantly by software patents. We wish to avoid the danger that redistributors of a free program will individually obtain patent licenses, in effect making the program proprietary. To prevent this, we have made it clear that any patent must be licensed for everyone's free use or not licensed at all.

The precise terms and conditions for copying, distribution and modification follow.

TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

    0. This License applies to any program or other work which contains a notice placed by the copyright holder saying it may be distributed under the terms of this Affero General Public License. The "Program", below, refers to any such program or work, and a "work based on the Program" means either the Program or any derivative work under copyright law: that is to say, a work containing the Program or a portion of it, either verbatim or with modifications and/or translated into another language. (Hereinafter, translation is included without limitation in the term "modification".) Each licensee is addressed as "you".

    Activities other than copying, distribution and modification are not covered by this License; they are outside its scope. The act of running the Program is not restricted, and the output from the Program is covered only if its contents constitute a work based on the Program (independent of having been made by running the Program). Whether that is true depends on what the Program does.
    1. You may copy and distribute verbatim copies of the Program's source code as you receive it, in any medium, provided that you conspicuously and appropriately publish on each copy an appropriate copyright notice and disclaimer of warranty; keep intact all the notices that refer to this License and to the absence of any warranty; and give any other recipients of the Program a copy of this License along with the Program.

    You may charge a fee for the physical act of transferring a copy, and you may at your option offer warranty protection in exchange for a fee.
    2. You may modify your copy or copies of the Program or any portion of it, thus forming a work based on the Program, and copy and distribute such modifications or work under the terms of Section 1 above, provided that you also meet all of these conditions:
        a) You must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.
        b) You must cause any work that you distribute or publish, that in whole or in part contains or is derived from the Program or any part thereof, to be licensed as a whole at no charge to all third parties under the terms of this License.
        c) If the modified program normally reads commands interactively when run, you must cause it, when started running for such interactive use in the most ordinary way, to print or display an announcement including an appropriate copyright notice and a notice that there is no warranty (or else, saying that you provide a warranty) and that users may redistribute the program under these conditions, and telling the user how to view a copy of this License. (Exception: if the Program itself is interactive but does not normally print such an announcement, your work based on the Program is not required to print an announcement.)
        d) If the Program as you received it is intended to interact with users through a computer network and if, in the version you received, any user interacting with the Program was given the opportunity to request transmission to that user of the Program's complete source code, you must not remove that facility from your modified version of the Program or work based on the Program, and must offer an equivalent opportunity for all users interacting with your Program through a computer network to request immediate transmission by HTTP of the complete source code of your modified version or other derivative work.

    These requirements apply to the modified work as a whole. If identifiable sections of that work are not derived from the Program, and can be reasonably considered independent and separate works in themselves, then this License, and its terms, do not apply to those sections when you distribute them as separate works. But when you distribute the same sections as part of a whole which is a work based on the Program, the distribution of the whole must be on the terms of this License, whose permissions for other licensees extend to the entire whole, and thus to each and every part regardless of who wrote it.

    Thus, it is not the intent of this section to claim rights or contest your rights to work written entirely by you; rather, the intent is to exercise the right to control the distribution of derivative or collective works based on the Program.

    In addition, mere aggregation of another work not based on the Program with the Program (or with a work based on the Program) on a volume of a storage or distribution medium does not bring the other work under the scope of this License.
    3. You may copy and distribute the Program (or a work based on it, under Section 2) in object code or executable form under the terms of Sections 1 and 2 above provided that you also do one of the following:
        a) Accompany it with the complete corresponding machine-readable source code, which must be distributed under the terms of Sections 1 and 2 above on a medium customarily used for software interchange; or,
        b) Accompany it with a written offer, valid for at least three years, to give any third party, for a charge no more than your cost of physically performing source distribution, a complete machine-readable copy of the corresponding source code, to be distributed under the terms of Sections 1 and 2 above on a medium customarily used for software interchange; or,
        c) Accompany it with the information you received as to the offer to distribute corresponding source code. (This alternative is allowed only for noncommercial distribution and only if you received the program in object code or executable form with such an offer, in accord with Subsection b above.)

    The source code for a work means the preferred form of the work for making modifications to it. For an executable work, complete source code means all the source code for all modules it contains, plus any associated interface definition files, plus the scripts used to control compilation and installation of the executable. However, as a special exception, the source code distributed need not include anything that is normally distributed (in either source or binary form) with the major components (compiler, kernel, and so on) of the operating system on which the executable runs, unless that component itself accompanies the executable.

    If distribution of executable or object code is made by offering access to copy from a designated place, then offering equivalent access to copy the source code from the same place counts as distribution of the source code, even though third parties are not compelled to copy the source along with the object code.
    4. You may not copy, modify, sublicense, or distribute the Program except as expressly provided under this License. Any attempt otherwise to copy, modify, sublicense or distribute the Program is void, and will automatically terminate your rights under this License. However, parties who have received copies, or rights, from you under this License will not have their licenses terminated so long as such parties remain in full compliance.
    5. You are not required to accept this License, since you have not signed it. However, nothing else grants you permission to modify or distribute the Program or its derivative works. These actions are prohibited by law if you do not accept this License. Therefore, by modifying or distributing the Program (or any work based on the Program), you indicate your acceptance of this License to do so, and all its terms and conditions for copying, distributing or modifying the Program or works based on it.
    6. Each time you redistribute the Program (or any work based on the Program), the recipient automatically receives a license from the original licensor to copy, distribute or modify the Program subject to these terms and conditions. You may not impose any further restrictions on the recipients' exercise of the rights granted herein. You are not responsible for enforcing compliance by third parties to this License.
    7. If, as a consequence of a court judgment or allegation of patent infringement or for any other reason (not limited to patent issues), conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License, they do not excuse you from the conditions of this License. If you cannot distribute so as to satisfy simultaneously your obligations under this License and any other pertinent obligations, then as a consequence you may not distribute the Program at all. For example, if a patent license would not permit royalty-free redistribution of the Program by all those who receive copies directly or indirectly through you, then the only way you could satisfy both it and this License would be to refrain entirely from distribution of the Program.

    If any portion of this section is held invalid or unenforceable under any particular circumstance, the balance of the section is intended to apply and the section as a whole is intended to apply in other circumstances.

    It is not the purpose of this section to induce you to infringe any patents or other property right claims or to contest validity of any such claims; this section has the sole purpose of protecting the integrity of the free software distribution system, which is implemented by public license practices. Many people have made generous contributions to the wide range of software distributed through that system in reliance on consistent application of that system; it is up to the author/donor to decide if he or she is willing to distribute software through any other system and a licensee cannot impose that choice.

    This section is intended to make thoroughly clear what is believed to be a consequence of the rest of this License.
    8. If the distribution and/or use of the Program is restricted in certain countries either by patents or by copyrighted interfaces, the original copyright holder who places the Program under this License may add an explicit geographical distribution limitation excluding those countries, so that distribution is permitted only in or among countries not thus excluded. In such case, this License incorporates the limitation as if written in the body of this License.
    9. Affero Inc. may publish revised and/or new versions of the Affero General Public License from time to time. Such new versions will be similar in spirit to the present version, but may differ in detail to address new problems or concerns.

    Each version is given a distinguishing version number. If the Program specifies a version number of this License which applies to it and "any later version", you have the option of following the terms and conditions either of that version or of any later version published by Affero, Inc. If the Program does not specify a version number of this License, you may choose any version ever published by Affero, Inc.

    You may also choose to redistribute modified versions of this program under any version of the Free Software Foundation's GNU General Public License version 3 or higher, so long as that version of the GNU GPL includes terms and conditions substantially equivalent to those of this license.
    10. If you wish to incorporate parts of the Program into other free programs whose distribution conditions are different, write to the author to ask for permission. For software which is copyrighted by Affero, Inc., write to us; we sometimes make exceptions for this. Our decision will be guided by the two goals of preserving the free status of all derivatives of our free software and of promoting the sharing and reuse of software generally.

    NO WARRANTY
    11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
    12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.)";

// https://stackoverflow.com/questions/13808714/accessing-individual-characters-in-a-file-inefficient-c
unsigned int countSequences(string fileName) {

	FILE *infile = fopen(fileName.c_str(), "r");

	char buffer[4096 * 1024];
	int size;
	unsigned int count = 0;

	int first = 0;
	int second = 1;

	for (int i = 0; i < size; i++) {
		if (buffer[i] == '>') {
			++count;
		}
	}

	while (0 < (size = fread(buffer, 1, sizeof(buffer), infile))) {
		for (int i = 0; i < size; i++)
			if (buffer[i] == '>')
				++count;
	}

	fclose(infile);

	return count;
}

double twoMeansAndDeviations(vector<double> &l, vector<double> &r, double f) {
	// Precondition
	if (r.size() != 4) {
		std::cerr << "twoMeansAndDeviations is expecting r of size 4, ";
		std::cerr << "but got " << r.size() << " elements.";
		std::cerr << std::endl;
		throw std::exception();
	}

	double m1 = r[0];
	double m2 = r[1];
	double s1 = r[2];
	double s2 = r[3];

	if (s1 <= 0.005) {
		s1 = 0.005;
	}
	if (s2 <= 0.005) {
		s2 = 0.005;
	}

	//	std::cout << "Initialization: " << m1 << " " << m2 << " " << s1 << " " << s2
	//			<< endl;

	vector<double> l1;
	vector<double> l2;
	int noChangeCounter = 0;
	double history1 = -1.0;
	double history2 = -1.0;
	for (int h = 0; h < 100; h++) {
		l1.clear();
		l2.clear();

		for (double id : l) {
			if (abs(m1 - id) / s1 < abs(m2 - id) / s2) {
				l1.push_back(id);
			} else {
				l2.push_back(id);
			}
		}

		if (l1.size() <= 1 || l2.size() <= 1) {
			// std::cout << "Converged. " << std::endl;
			break;
		}

		m1 = Util::calculateMean(l1);
		m2 = Util::calculateMean(l2);

		s1 = Util::calculateSTDev(l1, m1);
		s2 = Util::calculateSTDev(l2, m2);

		if (Util::isEqual(m1, history1) && Util::isEqual(m2, history2)) {
			noChangeCounter++;
			if (noChangeCounter == 3) {
				//				std::cout << "Stopping because there is no change for three iterations: ";
				//				std::cout << h << std::endl;
				break;
			}
		}

		history1 = m1;
		history2 = m2;

	}

	//	std::cout << m1 << " " << m2 << " " << s1 << " " << s2 << " "
	//			<< (double) l1.size() / l.size() << " "
	//			<< (double) l2.size() / l.size() << std::endl;

	double std1 = Util::calculateSTDev(l1, m1);
	double std2 = Util::calculateSTDev(l2, m2);

	double p1 = (double) l1.size() / l.size();
	double p2 = (double) l2.size() / l.size();

	double guess;
	if (p1 > p2) {
		guess = m1 - f * std1;
	} else {
		guess = m2 - f * std2;
	}
	return guess;
}

vector<double> twoMeans(vector<double> &l) {
	double m1 = *std::min_element(l.begin(), l.end());
	double m2 = *std::max_element(l.begin(), l.end());

	// std::cout << "Initialization: " << m1 << " " << m2 << endl;

	vector<double> l1;
	vector<double> l2;
	int noChangeCounter = 0;
	double history1 = -1.0;
	double history2 = -1.0;
	for (int h = 0; h < 100; h++) {
		l1.clear();
		l2.clear();

		for (double id : l) {
			if (abs(m1 - id) < abs(m2 - id)) {
				l1.push_back(id);
			} else {
				l2.push_back(id);
			}
		}

		if (l1.size() <= 1 || l2.size() <= 1) {
			// std::cout << "Converged. " << std::endl;
			break;
		}

		m1 = Util::calculateMean(l1);
		m2 = Util::calculateMean(l2);

		if (Util::isEqual(m1, history1) && Util::isEqual(m2, history2)) {
			noChangeCounter++;
			if (noChangeCounter == 3) {
				// std::cout << "Stopping because there is no change for three iterations: ";
				// std::cout << h << std::endl;
				break;
			}
		}

		history1 = m1;
		history2 = m2;
	}

	//	std::cout << m1 << " " << m2 << " " << (double) l1.size() / l.size() << " "
	//			<< (double) l2.size() / l.size() << std::endl;

	// Fill the results
	vector<double> r(4, 0.0);
	r[0] = m1;
	r[1] = m2;
	r[2] = Util::calculateSTDev(l1, m1);
	r[3] = Util::calculateSTDev(l2, m2);
	return r;
}

/**
 * Estimate the bandwidth based on one block
 */
template<class V>
vector<double> guessBandwidthOneTime(Block *block, int cores,
		IdentityCalculator<V> &identity) {
	int size = block->size();
	auto tup = identity.unpackBlock(block, cores);
	V **kHistList = std::get<0>(tup);
	uint64_t **monoHistList = std::get<1>(tup);
	std::string **infoList = std::get<2>(tup);
	int *lenList = std::get<3>(tup);

	Matrix a = identity.score(kHistList, monoHistList, size, cores, lenList);

	int r = a.getNumRow();

	// Collect the 4th closest sequence to every sequence,
	// assuming the minimum cluster size to be 5
	int minSize = 5;
	vector<double> v(r * minSize, -1.0);

#pragma omp parallel for schedule(static) num_threads(cores)
	for (int i = 0; i < r; i++) {
		vector<double> row = a.getRow(i);
		std::nth_element(row.begin(), row.begin() + (minSize - 1), row.end(),
				std::greater { });
		for (int h = 0; h < minSize; h++) {
			v[i * minSize + h] = row[h];
		}
	}

	// Remove r ones representing the score of a sequence vs. itself
	vector<double> vNoOnes; // It may contain ones but not due to a sequence vs. itself
	int h = 0;
	for (double d : v) {
		if (h >= r) {
			vNoOnes.push_back(d);
		} else {
			if (d == 1.0) {
				h++;
			} else {
				vNoOnes.push_back(d);
			}
		}
	}

	return vNoOnes;
}

template<class V>
double guessBandwidthHelper(string dbFile, int cores, DataGenerator *g) {
	// Count sequences
	unsigned int seqCount = countSequences(dbFile);

	double ratio = (double) Parameters::getMsBandwidthBlock() / seqCount;
	double sigmas = 0.0; // < 5%
	if (ratio > 0.05 && ratio <= 0.25) {
		// Between 5% and 25%
		sigmas = 1.0;
	} else if (ratio > 0.25) {
		// More than 25%
		sigmas = 2.0;
	}

	// Train Identity
	IdentityCalculator<V> identity(g, cores,
			Parameters::getMsBandwidthThreshold(), false /*canSkip*/,
			false /*canRelax*/);

	// Estimate on a small number of blocks
	vector<double> guessList;

	FastaReader reader(dbFile, Parameters::getMsBandwidthBlock());

	if (ratio < 1.0) {
		for (int i = 0;
				i < Parameters::getMsBandwidthIterations()
						&& reader.isStillReading(); i++) {
			Block *block = reader.read();

			// Skip the last block if it is small
			if (i > 0 && block->size() < Parameters::getMsBandwidthBlock()) {
				break;
			}
			vector<double> idList = guessBandwidthOneTime<V>(block, cores,
					identity);

			vector<double> r = twoMeans(idList);
			double guess = twoMeansAndDeviations(idList, r, sigmas);
			cout << "Estimated threshold " << i << ": " << guess << endl;
			guessList.push_back(guess);
		}
	} else {
		Block *block = reader.read();
		vector<double> idList = guessBandwidthOneTime<V>(block, cores,
				identity);

		double m = Util::calculateMean(idList);
		double std = Util::calculateSTDev(idList, m);
		double min = *std::min_element(idList.begin(), idList.end());
		//		cout << "Mean = " << m << endl;
		//		cout << "STD = " << std << endl;
		//		cout << "Min = " << min << endl;
		double guess = m - 3.0 * std;
		if (guess < min) {
			guess = min;
		}
		guessList.push_back(guess);
	}

	// The threshold is the median (the average if only two trials)
	double threshold = -1.0;
	if (guessList.size() > 2) {
		int medianIndex = guessList.size() / 2;
		std::nth_element(guessList.begin(), guessList.begin() + medianIndex,
				guessList.end(), std::greater { });
		threshold = guessList.at(medianIndex);
	} else if (guessList.size() == 2) {
		threshold = (guessList.at(0) + guessList.at(1)) / 2.0;
	} else if (guessList.size() == 1) {
		threshold = guessList.at(0);
	}

	// Post-condition
	if (threshold == -1.0) {
		std::cerr << "Something wrong with determining the threshold.";
		std::cerr << std::endl;
		throw std::exception();
	}

	// For testing only
	//	std::cout << "============================================" << std::endl;
	//	for (double g : guessList) {
	//		std::cout << g << std::endl;
	//	}
	// End testing
	// std::cout << "\tFinal threshold: " << threshold << std::endl;

	return threshold;
}

double guessBandwidth(string dbFile, int cores) {
	std::cout << "Estimating the threshold ..." << std::endl;
	SynDataGenerator g(dbFile, Parameters::getMsBandwidthThreshold(), cores);
	int64_t maxLength = g.getMaxLength();
	double threshold = -1.0;
	// Determine histogram data type
	if (maxLength <= std::numeric_limits<int8_t>::max()) {
		threshold = guessBandwidthHelper<int8_t>(dbFile, cores, &g);
	} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
		threshold = guessBandwidthHelper<int16_t>(dbFile, cores, &g);
	} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
		threshold = guessBandwidthHelper<int32_t>(dbFile, cores, &g);
	} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
		threshold = guessBandwidthHelper<int64_t>(dbFile, cores, &g);
	} else {
		std::cout << "Warning: Overflow is possible however unlikely.";
		std::cout << std::endl;
		threshold = guessBandwidthHelper<int64_t>(dbFile, cores, &g);
	}

	if (threshold < 0 || threshold > 1.0) {
		std::cerr << "A threshold must be between 0 and 1." << std::endl;
		std::cerr << "Could not determine the threshold." << std::endl;
		throw std::exception();
	}
	return threshold;
}

template<class V>
void start(string dbFile, int blockSize, int vBlockSize, int passNum,
		DataGenerator *g, int threadNum, double t, bool canAssignAll,
		string outFile, bool canEvaluate,
		bool msCanRelax /*Relax the final assignment to include extended members*/) {

	bool canSkip = true;

	if (t > 0.99) {
		bool idCanRelax = false;
		double idThreshold = t > 0.99 ? 0.99 : t;
		IdentityCalculator1<V> identity(g, threadNum, idThreshold, canSkip,
				idCanRelax);
		MeanShiftLarge<V> ms(dbFile, blockSize, vBlockSize, passNum, identity,
				threadNum, canAssignAll, outFile, t, canEvaluate, msCanRelax);
	} else {
		bool idCanRelax = true;
		IdentityCalculator<V> identity(g, threadNum, t, canSkip, idCanRelax);
		MeanShiftLarge<V> ms(dbFile, blockSize, vBlockSize, passNum, identity,
				threadNum, canAssignAll, outFile, t, canEvaluate, msCanRelax);
	}
}

int main(int argc, char *argv[]) {
	std::cout << std::endl;
	std::cout << "MeShClust v3.0 is developed by Hani Z. Girgis, PhD."
			<< std::endl;
	std::cout << std::endl;
	std::cout
			<< "This program clusters DNA sequences using identity scores obtained without alignment."
			<< std::endl;
	std::cout << std::endl;
	std::cout << "Copyright (C) 2021-2022 Hani Z. Girgis, PhD" << std::endl;
	std::cout << std::endl;
	std::cout << "Academic use: Affero General Public License version 1."
			<< std::endl;
	std::cout << std::endl;
	std::cout
			<< "Any restrictions to use for profit or non-academics: Alternative commercial license is required."
			<< std::endl;
	std::cout << std::endl;
	std::cout
			<< "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;"
			<< std::endl;
	std::cout
			<< "without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE."
			<< std::endl;
	std::cout << std::endl;
	std::cout
			<< "Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information."
			<< std::endl;
	std::cout << std::endl;
	std::cout << "Please cite the following papers: " << std::endl;

	std::cout << "\t"
			<< "1. Identity: Rapid alignment-free prediction of sequence alignment identity scores using"
			<< std::endl;
	std::cout << "\t"
			<< "self-supervised general linear models. Hani Z. Girgis, Benjamin T. James, and Brian B."
			<< std::endl;
	std::cout << "\t" << "Luczak. NAR GAB, 3(1):lqab001, 2021." << std::endl;

	std::cout << "\t"
			<< "2. MeShClust: an intelligent tool for clustering DNA sequences. Benjamin T. James,"
			<< std::endl;
	std::cout << "\t"
			<< "Brian B. Luczak, and Hani Z. Girgis. Nucleic Acids Res, 46(14):e83, 2018."
			<< std::endl;

	std::cout << "\t"
			<< "3. MeShClust v3.0: High-quality clustering of DNA sequences using the mean shift algorithm"
			<< std::endl
			<< "\tand alignment-free identity scores. Hani Z. Girgis. "
			<< "A great journal. " << "2022." << std::endl;

	std::cout << std::endl;

	if (argc == 1 || (argc == 2 && argv[1][1] == 'h')) {
		std::cout << "List of parameters:" << std::endl;
		// Required parameters
		std::cout << "\t-d: Required. Database file in FASTA format."
				<< std::endl;
		std::cout
				<< "\t-o: Required. Output file. Each line has 4 tab-separated fields: cluster number, sequence header,"
				<< std::endl;
		std::cout
				<< "\t    identity score with the cluster center, C/M/E/O. C/M/E/O stand for center, member, extended"
				<< std::endl;
		std::cout
				<< "\t    member (threshold - regression error), outside (less than threshold). The O mark should be seen"
				<< std::endl;
		std::cout
				<< "\t    when the -a y is used."
				<< std::endl;

		// Optional parameters
		std::cout
				<< "\t-t: Optional. Threshold identity score (between 0 & 0.99) for determining cluster membership."
				<< std::endl;

		std::cout
				<< "\t-a: Optional. Assign every sequence to a cluster regardless of the threshold -- y or n"
				<< std::endl;
		std::cout
				<< "\t    (default: n). If no, a sequence that is not within the threshold score of any"
				<< std::endl;
		std::cout
				<< "\t    cluster will comprise its own cluster. If yes, the assignment step may take long time on large sets."
				<< std::endl;
		std::cout
				<< "\t    It would not take additional time if the evaluation option and this option are enabled together."
				<< std::endl;

		std::cout
				<< "\t-c: Optional. Number of cores or hyperthreads. For the search mode, set this parameter to the"
				<< std::endl;
		std::cout
				<< "\t    number of cores not hyperthreads. For example, suppose your computer has 4 cores, each of"
				<< std::endl;
		std::cout
				<< "\t    which supports 2 hyperthreads. Set this parameter to 4 if you are using the search mode or"
				<< std::endl;
		std::cout
				<< "\t    to 8 if you are using the all-versus-all mode. By default, all hyperthreads are used."
				<< std::endl;

		std::cout
				<< "\t-r: Optional. Automatically relax the threshold according to the predictor error -- y or n"
				<< std::endl;
		std::cout
				<< "\t    (default: y). This option affects the final assignment step only."
				<< std::endl;

		std::cout
				<< "\t-e: Optional. Evaluate cluster quality. May take long time on large data sets -- y or n"
				<< std::endl;
		std::cout
				<< "\t    (default: n). It would not take additional time if this option and the assign-all "
				<< std::endl;
		std::cout << "\t    option are enabled together." << std::endl;

		std::cout
				<< "\t-b: Optional. The batch size for all vs. all (default: 25,000; maximum: 46,340)."
				<< std::endl;
		std::cout << "\t    Increasing this number will slow the program."
				<< std::endl;

		std::cout
				<< "\t-v: Optional. The batch size of sequences to be read (default: 100,000)."
				<< std::endl;
		std::cout
				<< "\t    It is recommended to be 2-4 times the all-vs-all batch adjusted by parameter -b."
				<< std::endl;
		std::cout << "\t    Increasing this number will require more memory."
				<< std::endl;

		std::cout << "\t-p: Optional. The number of data passes (default: 10)."
				<< std::endl;
		std::cout
				<< "\t    It applies to the scaled-up version -- not to the original algorithm."
				<< std::endl;

		std::cout
				<< "\t-l: Optional. Print academic license (Affero General Public License version 1) and exit -- y"
				<< std::endl;
		std::cout << "\t    (yes) or n (no)." << std::endl;
		std::cout << "\t-h: Optional. Print this help message." << std::endl;

		std::cout << std::endl;

		std::cout << "Examples: " << std::endl;
		std::cout
				<< "\t1. To cluster sequences with a minimum identity score of 0.8"
				<< std::endl;
		std::cout << "\t\tmeshclust -d input.fa -o output.txt -t 0.8"
				<< std::endl;
		std::cout << std::endl;

		std::cout
				<< "\t2. To cluster sequences with an estimated minimum identity score"
				<< std::endl;
		std::cout << "\t\tmeshclust -d input.fa -o output.txt" << std::endl;
		std::cout << std::endl;

		std::cout
				<< "\t3. To cluster sequences with a minimum identity score of 0.8 and evaluate"
				<< std::endl;
		std::cout << "\t\tmeshclust -d input.fa -o output.txt -t 0.8 -e y"
				<< std::endl;
		std::cout << std::endl;

		std::cout
				<< "\t4. To cluster sequences with a minimum identity score of 0.8 and assign every"
				<< std::endl;
		std::cout
				<< "\t   sequence to a cluster even if its identity score with a cluster center is"
				<< std::endl;
		std::cout
				<< "\t   less than the minimum score (useful when your data are noise free)"
				<< std::endl;
		std::cout << "\t\tmeshclust -d input.fa -o output.txt -t 0.8 -a y"
				<< std::endl;
		std::cout << std::endl;

		std::cout
				<< "\t5. To cluster sequences with a minimum identity score of 0.8 and assign each sequence"
				<< std::endl;
		std::cout
				<< "\t   according to the minimum score without relaxing by the regression model error"
				<< std::endl;
		std::cout << "\t\tmeshclust -d input.fa -o output.txt -t 0.8 -r n"
				<< std::endl;
		std::cout << std::endl;

		std::cout
				<< "\t6. To cluster sequences with a minimum identity score of 0.8 using an all-vs-all"
				<< std::endl;
		std::cout << "\t   block size of 1000 and a reading block size of 4000"
				<< std::endl;

		std::cout
				<< "\t\tmeshclust -d input.fa -o output.txt -t 0.8 -b 1000 -v 4000"
				<< std::endl;
		std::cout << std::endl;

		std::cout
				<< "\t7. To cluster sequences with a minimum identity score of 0.8 and specify the number"
				<< std::endl;
		std::cout
				<< "\t   of data passes (useful when the algorithm did not converge, i.e., cluster count"
				<< std::endl;
		std::cout << "\t   kept changing from iteration to iteration)"
				<< std::endl;
		std::cout << "\t\tmeshclust -d input.fa -o output.txt -t 0.8 -p 100"
				<< std::endl;
		std::cout << std::endl;

		std::cout << "\t8. To print the academic license" << std::endl;
		std::cout << "\t\tmeshclust -l y" << std::endl;
		std::cout << std::endl;

		exit(0);
	}

	std::string dbFile;
	std::string qryFile;
	std::string outFile;

	char relax = 'y';
	char evaluate = 'n';
	char license = 'n';
	char all = 'n';
	int cores = std::thread::hardware_concurrency();
	double threshold = 0.0;
	bool isThresholdProvided = false;
	int blockSize = 0;
	bool isBlockSizeProvided = false;
	int vBlockSize = 0;
	bool isVBlockSizeProvided = false;
	int passNum = 0;
	bool isPassNumProvided = false;

	for (int i = 1; i < argc; i += 2) {
		switch (argv[i][1]) {
		case 'd': {
			dbFile = std::string(argv[i + 1]);
		}
			break;

		case 'o': {
			outFile = std::string(argv[i + 1]);
		}
			break;

		case 'c': {
			cores = atoi(argv[i + 1]);
		}
			break;

		case 't': {
			threshold = atof(argv[i + 1]);
			isThresholdProvided = true;
		}
			break;

		case 'r': {
			relax = argv[i + 1][0];
		}
			break;

		case 'e': {
			evaluate = argv[i + 1][0];
		}
			break;

		case 'l': {
			license = argv[i + 1][0];
		}
			break;

		case 'a': {
			all = argv[i + 1][0];
		}
			break;

		case 'b': {
			blockSize = atoi(argv[i + 1]);
			isBlockSizeProvided = true;
		}
			break;

		case 'v': {
			vBlockSize = atoi(argv[i + 1]);
			isVBlockSizeProvided = true;
		}
			break;

		case 'p': {
			passNum = atoi(argv[i + 1]);
			isPassNumProvided = true;
		}
			break;

		default: {
			std::cerr << argv[i][1]
					<< " is invalid option. Rerun with -h to see the help message.";
			std::cerr << std::endl;
			std::cerr << std::endl;
			exit(1);
		}
		}
	}

	if (license == 'y') {
		std::cout << agplv1 << std::endl;
		exit(0);
	} else if (license == 'n') {
		// Do nothing
	} else {
		std::cerr
				<< "Error: If you would like to print the academic license use -l y, otherwise -l n.";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	if (all != 'y' && all != 'n') {
		std::cerr
				<< "Error: If you would like to print the scores of all pairs use -a y, otherwise -a n.";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	// Make sure that the required parameters have been provided
	if (dbFile.empty()) {
		std::cerr << "Error: Please provide a database file in FASTA format.";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	if (outFile.empty()) {
		std::cerr << "Error: Please provide an output file.";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	if (isThresholdProvided && (threshold < 0.0 || threshold > 1.0)) {
		std::cerr << "Error: Please provide a threshold between 0.00 and 1.00";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	if (relax != 'y' && relax != 'n') {
		std::cerr
				<< "Error: Please provide a valid answer -- y (relax threshold) or n (do not relax threshold).";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	if (evaluate != 'y' && evaluate != 'n') {
		std::cerr
				<< "Error: Please provide a valid answer -- y (evaluate cluster quality) or n (do not evaluate).";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	if (isBlockSizeProvided
			&& (blockSize < 1000 || blockSize > Parameters::getMsMaxMatrixSize())) {
		std::cerr
				<< "Error: Please provide a batch (all vs all) size between 1,000 and 46,340.";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	if (isVBlockSizeProvided) {
		if ((isBlockSizeProvided && vBlockSize < blockSize)
				|| (!isBlockSizeProvided
						&& vBlockSize < Parameters::getMsBlock())) {
			std::cerr
					<< "Error: Please provide a read-sequences batch size >= the all-vs-all batch size.";
			std::cerr << std::endl;
			std::cerr << "\tRerun with -h to see the help message.";
			std::cerr << std::endl;
			std::cerr << std::endl;
			exit(1);
		}
	}

	if (isPassNumProvided && (passNum < 2 || passNum > 10)) {
		std::cerr
				<< "Error: Please provide a number for data passes between 2 and 10.";
		std::cerr << std::endl;
		std::cerr << "\tRerun with -h to see the help message.";
		std::cerr << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	// Print parameter values
	std::cout << "Database file: " << dbFile << std::endl;
	std::cout << "Output file: " << outFile << std::endl;
	std::cout << "Cores: " << cores << std::endl;

	if (isThresholdProvided) {
		std::cout << "Provided threshold: " << threshold << std::endl;
	} else if (threshold == 0.0) {
		threshold = guessBandwidth(dbFile, cores);
		std::cout << "Calculated threshold: " << threshold << std::endl;
	}

	if (blockSize == 0) {
		blockSize = Parameters::getMsBlock();
	}
	std::cout << "Block size for all vs. all: " << blockSize << std::endl;

	if (vBlockSize == 0) {
		vBlockSize = Parameters::getMsVBlock();
	}
	std::cout << "Block size for reading sequences: " << vBlockSize
			<< std::endl;

	if (passNum == 0) {
		passNum = Parameters::getMsPassNum();
	}
	std::cout << "Number of data passes: " << passNum << std::endl;

	bool canAssignAll = (all == 'y') ? true : false;
	std::cout << "Can assign all: " << (canAssignAll ? "Yes" : "No")
			<< std::endl;
	std::cout << std::endl << std::endl;

	SynDataGenerator g(dbFile, threshold, cores);
	int64_t maxLength = g.getMaxLength();

	bool canEvaluate = (evaluate == 'y') ? true : false;
	bool canRelax = (relax == 'y') ? true : false;

	// Determine histogram data type
	if (maxLength <= std::numeric_limits<int8_t>::max()) {
		start<int8_t>(dbFile, blockSize, vBlockSize, passNum, &g, cores,
				threshold, canAssignAll, outFile, canEvaluate, canRelax);
	} else if (maxLength <= std::numeric_limits<int16_t>::max()) {
		start<int16_t>(dbFile, blockSize, vBlockSize, passNum, &g, cores,
				threshold, canAssignAll, outFile, canEvaluate, canRelax);
	} else if (maxLength <= std::numeric_limits<int32_t>::max()) {
		start<int32_t>(dbFile, blockSize, vBlockSize, passNum, &g, cores,
				threshold, canAssignAll, outFile, canEvaluate, canRelax);
	} else if (maxLength <= std::numeric_limits<int64_t>::max()) {
		start<int64_t>(dbFile, blockSize, vBlockSize, passNum, &g, cores,
				threshold, canAssignAll, outFile, canEvaluate, canRelax);
	} else {
		std::cout << "ReaderAlignerCoordinator warning: ";
		std::cout << "Overflow is possible however unlikely.";
		std::cout << std::endl;
		std::cout << "A histogram entry is 64 bits." << std::endl;
		start<int64_t>(dbFile, blockSize, vBlockSize, passNum, &g, cores,
				threshold, canAssignAll, outFile, canEvaluate, canRelax);
	}

	std::cout << "Finished." << std::endl;
	std::cout << std::endl;
	std::cout
			<< "Thanks for using MeShClust v3.0. Please post any questions or problems on GitHub: "
			<< std::endl;
	std::cout
			<< "https://github.com/BioinformaticsToolsmith/Identity or email Dr. Hani Z. Girgis."
			<< std::endl;
	std::cout << std::endl;
	return 0;
}
