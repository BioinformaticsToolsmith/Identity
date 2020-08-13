# Identity

Identity 1.0 is developed by Hani Z. Girgis, PhD.

This program calculates DNA sequence identity scores rapidly without alignment.

Copyright (C) 2020 Hani Z. Girgis, PhD

Academic use: Affero General Public License version 1.

Any restrictions to use for profit or non-academics: Alternative commercial license is required.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.

Please cite the following paper: 

	Identity: Rapid alignment-free prediction of sequence alignment identity scores using
	self-supervised general linear models. Hani Z. Girgis, Benjamin T. James, and Brian B.
	Luczak. NAR GAB, 2020.

Requirments:

	GNU g++ 7.5.0 or later
	You may change the compiler in the CMakeLists.txt: set(CMAKE_CXX_COMPILER g++-7)

To Compile:

	mkdir bin
	cd bin
	cmake ..
	make

To Test:

	cd test
	../bin/identity -d keratin_small.fasta -q keratin_query.fasta -o output.txt -t 0.7

List of parameters:

	-d: Required. Database file in FASTA format.
	
	-o: Required. Output file. Each line has 3 tab-separated fields (>header1    >header2    score).
	
	-t: Required. Identity score threshold (between 0 & 0.99), below which pairs are not reported.
	
	-q: Optional. Query file in FASTA format. If no query(s) is provided, all versus all is
	    performed on the database file.
	    
	-c: Optional. Number of cores or hyperthreads. For the search mode, set this parameter to the
	    number of cores not hyperthreads. For example, suppose your computer has 4 cores, each of
	    which supports 2 hyperthreads. Set this parameter to 4 if you are using the search mode or
	    to 8 if you are using the all-versus-all mode. By default, all hyperthreads are used.
	    
	-r: Optional. Automatically relax the threshold according to the predictor error -- y (yes) or
	    n (no). By default, it is enabled except if the threshold is 0.9 or higher.
	    
	-l: Optional. Print academic license (Affero General Public License version 1) and exit -- y
	    (yes) or n (no).
	    
	-h: Optional. Print this help message.

Examples: 

	1. To perform database search with a minimum identity score of 0.7
		identity -d databas.fasta -q query.fasta -o output.txt -t 0.7

	2. To perform database search with a minimum identity score of 0.7 using 10 threads
		identity -d databas.fasta -q query.fasta -o output.txt -t 0.7 -c 10

	3. To perform all versus all with a minimum identity score of 0.8
		identity -d databas.fasta -o output.txt -t 0.8

	4. To perform all versus all with a minimum identity score of 0.8 with strict threshold
		identity -d databas.fasta -o output.txt -t 0.8 -r n

	5. To print the academic lincense
		identity -l y
		
Hierarchical clustering:

	To utilize the all-versus-all identity scores in hierarchical clustering see the README file 
	under the py directory.

