# Identity & MeShClust

Identity 1.1 is developed by Hani Z. Girgis, PhD.

This program calculates DNA sequence identity scores rapidly without alignment.

Copyright (C) 2020–2022 Hani Z. Girgis, PhD

Academic use: Affero General Public License version 1.

Any restrictions to use for profit or non-academics: Alternative commercial license is required.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.

Please cite the following paper: 

	Identity: Rapid alignment-free prediction of sequence alignment identity scores using 
	self-supervised general linear models (2021). Hani Z. Girgis, Benjamin T. James, and 
	Brian B. Luczak. NAR Genom Bioinform, 13(1), lqab001.

Requirments:

	GNU g++ 7.5.0 or later

To Compile:

	mkdir bin
	cd bin
	(If your default compiler meets the version requiremet) 
	cmake ..
	(Or if you would like to specify a different compiler that meets the requirement)
 	cmake .. -DCMAKE_CXX_COMPILER=your_compiler_name_for_example_g++-7
	make

# Identity

To Test:

	cd test
	../bin/identity -d keratin_small.fasta -q keratin_query.fasta -o output.txt -t 0.7

List of parameters:

	-d: Required. Database file in FASTA format.
	
	-o: Required. Output file. Each line has 3 tab-separated fields (>header1    >header2    score).
	
	-t: Required. Identity score threshold (between 0 & 0.99), below which pairs are not reported.
	
	-a: Optional. Report identity scores for all pairs including those below the threshold -- y
	    (yes) or n (no). If yes, it may take long time on large datasets due to writing to a file.
	    This option should be used if you desire constructing a phylogenetic tree. You may use the
	    accompanying Python program to covert Identity's output to a Phylip distance matrix.
	
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

	5. To perform all versus all with a minimum identity score of 0.8 and report all pairs
		identity -d databas.fasta -o output.txt -t 0.8 -a y

	6. To print the academic lincense
		identity -l y
		
Phylogenetic trees:

	To produce an all-versus-all distance matrix in Phylip format use makePhylipMatrix.py under the
	py directory. 

# MeShClust

To Test:

	cd test
	../bin/meshclust -d 97_shuffled.fa -o output.txt -t 0.97
	
	MeShClust should produce 100 clusters.

List of parameters:

	-d: Required. Database file in FASTA format.
	-o: Required. Output file. Each line has 4 tab-separated fields: cluster number, sequence header,
	    identity score with the cluster center, C/M/E/O. C/M/E/O stand for center, member, extended
	    member (threshold - regression error), outside (less than threshold). The O mark should be seen
	    when the -a y is used.
	-t: Optional. Threshold identity score (between 0 & 0.99) for determining cluster membership.
	-a: Optional. Assign every sequence to a cluster regardless of the threshold -- y or n
	    (default: n). If no, a sequence that is not within the threshold score of any
	    cluster will comprise its own cluster. If yes, the assignment step may take long time on large sets.
	    It would not take additional time if the evaluation option and this option are enabled together.
	-c: Optional. Number of cores or hyperthreads. For the search mode, set this parameter to the
	    number of cores not hyperthreads. For example, suppose your computer has 4 cores, each of
	    which supports 2 hyperthreads. Set this parameter to 4 if you are using the search mode or
	    to 8 if you are using the all-versus-all mode. By default, all hyperthreads are used.
	-r: Optional. Automatically relax the threshold according to the predictor error -- y or n
	    (default: y). This option affects the final assignment step only.
	-e: Optional. Evaluate cluster quality. May take long time on large data sets -- y or n
	    (default: n). It would not take additional time if this option and the assign-all 
	    option are enabled together.
	-b: Optional. The batch size for all vs. all (default: 25,000; maximum: 46,340).
	    Increasing this number will slow the program.
	-v: Optional. The batch size of sequences to be read (default: 100,000).
	    It is recommended to be 2-4 times the all-vs-all batch adjusted by parameter -b.
	    Increasing this number will require more memory.
	-p: Optional. The number of data passes (default: 10).
	    It applies to the scaled-up version -- not to the original algorithm.
	-l: Optional. Print academic license (Affero General Public License version 1) and exit -- y
	    (yes) or n (no).
	-h: Optional. Print this help message.

Output format:

	Each line of the output file has four fields.
	Field 1 is the cluster identifier.
	Field 2 is the sequence identifier.
	Field 3 is the identity score between a member sequence and the sequence representing the center of a cluster.
	Field 4 is the status of a sequence. One of four letters (C, M, E, O) may appear in this field. C stands for 
		the center of the cluster. M represents a member whose identity score with the center is less than or
		equal to the threshold. E stands for extended member whose identity score with the center is less 
		than or equal to the relaxed threshold (threshold – regression error). O stands for outside. The O 
		letter may appear when the assign-all option is enabled (-a y); in this case a sequence is assigned to 
		the closest center regardless of the threshold or the relaxed threshold.
		
Examples: 

	1. To cluster sequences with a minimum identity score of 0.8
		meshclust -d input.fa -o output.txt -t 0.8

	2. To cluster sequences with an estimated minimum identity score
		meshclust -d input.fa -o output.txt

	3. To cluster sequences with a minimum identity score of 0.8 and evaluate
		meshclust -d input.fa -o output.txt -t 0.8 -e y

	4. To cluster sequences with a minimum identity score of 0.8 and assign every
	   sequence to a cluster even if its identity score with a cluster center is
	   less than the minimum score (useful when your data are noise free)
		meshclust -d input.fa -o output.txt -t 0.8 -a y

	5. To cluster sequences with a minimum identity score of 0.8 and assign each sequence
	   according to the minimum score without relaxing by the regression model error
		meshclust -d input.fa -o output.txt -t 0.8 -r n

	6. To cluster sequences with a minimum identity score of 0.8 using an all-vs-all
	   block size of 1000 and a reading block size of 4000
		meshclust -d input.fa -o output.txt -t 0.8 -b 1000 -v 4000

	7. To cluster sequences with a minimum identity score of 0.8 and specify the number
	   of data passes (useful when the algorithm did not converge, i.e., cluster count
	   kept changing from iteration to iteration)
		meshclust -d input.fa -o output.txt -t 0.8 -p 100

	8. To print the academic license
		meshclust -l y


 
		 
