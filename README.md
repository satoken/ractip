RactIP for predicting RNA-RNA interaction using integer programming
===================================================================

Requirements
------------

* [Vienna RNA package](http://www.tbi.univie.ac.at/~ivo/RNA/) (>= 1.8)
* [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
  or [Gurobi Optimizer](http://www.gurobi.com/) (>=2.0)
  or [ILOG CPLEX](http://www.ibm.com/software/products/ibmilogcple/) (>=12.0)

Install
-------

For GLPK,

	./configure --with-vienna-rna=/path/to/vienna-rna --with-glpk

For Gurobi,

	./configure --with-vienna-rna=/path/to/vienna-rna --with-gurobi=/path/to/gurobi

For CPLEX,

	./configure --with-vienna-rna=/path/to/vienna-rna --with-cplex \
	            --with-cplex-include=/path/to/cplex/include \
		        --with-cplex-lib=/path/to/cplex/lib

Then,

	make
	make install

Usage
-----

RactIP can take two FASTA formatted RNA sequences as input, the
predict their joint secondary structures.

	Usage: ractip [OPTIONS]... [FILES]...

	-h, --help                 Print help and exit
        --full-help            Print help, including hidden options, and exit
	-V, --version              Print version and exit
	-a, --alpha=FLOAT          weight for hybridization  (default=`0.6')
	-b, --beta=FLOAT           weight for accessibility  (default=`0.0')
	-t, --fold-th=FLOAT        Threshold for base-pairing probabilities
                               (default=`0.8')
	-u, --hybridize-th=FLOAT   Threshold for hybridazation probabilities
                               (default=`0.3')
	-s, --acc-th=FLOAT         Threshold for accessible probabilities
                               (default=`0.005')
	-e, --show-energy          calculate the free energy of the predicted joint
                               structure  (default=off)

	% ractip DIS.fa DIS.fa
	>DIS
	CUCGGCUUGCUGAGGUGCACACAGCAAGAGGCGAG
	((((.(((((((..[[[[[[.)))))))...))))
	>DIS
	CUCGGCUUGCUGAGGUGCACACAGCAAGAGGCGAG
	((((.(((((((..]]]]]].)))))))...))))

The parenthesis '()' and the brackets '[]' indicate the predicted
internal base-pairs and external base-pairs (interactions),
respectively. 


References
----------

* Kato, Y., Sato, K., Hamada, M., Watanabe, Y., Asai, K., Akutsu, T.:
  RactIP: fast and accurate prediction of RNA-RNA interaction using
  integer programming. *Bioinformatics*, 26(18):i460-i466, 2010.

