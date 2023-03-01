RactIP for predicting RNA-RNA interaction using integer programming
===================================================================

Requirements
------------
* C++17 compatible compiler (tested on Apple clang version 14.0.0 and GCC version 10.2.1)
* cmake (>= 3.8)
* pkg-config
* [Vienna RNA package](https://www.tbi.univie.ac.at/RNA/) (>= 2.2.0)
* one of these MIP solvers
    * [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
    * [Gurobi Optimizer](http://www.gurobi.com/) (>=8.0)
    * [ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) (>=12.0)
    * [SCIP](https://scipopt.org/) (>= 8.0.3)
    * [HiGHS](https://highs.dev/) (>= 1.5.0)

Install
-------

For GLPK,

	export PKG_CONFIG_PATH=/path/to/viennarna/lib/pkgconfig:$PKG_CONFIG_PATH
	mkdir build && cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..  # configure
	cmake --build . # build
	cmake --install . # install (optional)

For Gurobi, add ``-DENABLE_GUROBI`` to the configure step:

	cmake -DENABLE_GUROBI -DCMAKE_BUILD_TYPE=Release ..  # configure

For CPLEX, add ``-DENABLE_CPLEX`` to the configure step:

	cmake -DENABLE_CPLEX -DCMAKE_BUILD_TYPE=Release ..  # configure

For SCIP, add ``-DENABLE_SCIP`` to the configure step:

    cmake -DENABLE_SCIP -DCMAKE_BUILD_TYPE=Release ..  # configure

For HiGHS, add ``-DENABLE_HIGHS`` to the configure step:

    cmake -DENABLE_HIGHS -DCMAKE_BUILD_TYPE=Release ..  # configure

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

### Run with Docker

    docker build . -t ractip
    docker run -it --rm -v $(pwd):$(pwd) -w $(pwd) ractip ractip DIS.fa DIS.fa

### Run with the web server

The RactIP web server is available at http://ws.sato-lab.org/rtips/ractip/.


References
----------

* Kato, Y., Sato, K., Hamada, M., Watanabe, Y., Asai, K., Akutsu, T.:
  RactIP: fast and accurate prediction of RNA-RNA interaction using
  integer programming. *Bioinformatics*, 26(18):i460-i466, 2010.
* Kato, Y., Mori, T., Sato, K., Maegawa, S., Hosokawa, H., Akutsu, T.:
  An accessibility-incorporated method for accurate prediction of RNA–RNA 
  interactions from sequence data, *Bioinformatics*, 33(2):202-209, 2017.
