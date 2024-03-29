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
	cmake -DCMAKE_BUILD_TYPE=Release -S . -B build  # configure
	cmake --build build  # build
	sudo cmake --install build  # install (optional)

If GLPK solver has been installed to the directory other than ``/usr`` or ``/usr/local``, specify the option ``-DGLPK_ROOT_DIR=/path/to/glpk`` in the configure step.

For Gurobi, add ``-DENABLE_GUROBI`` to the configure step:

	cmake -DENABLE_GUROBI=true -DGUROBI_DIR=/path/to/gurobi -DCMAKE_BUILD_TYPE=Release -S . -B build  # configure

For CPLEX, add ``-DENABLE_CPLEX`` to the configure step:

	cmake -DENABLE_CPLEX=true -DCPLEX_ROOT_DIR=/path/to/cplex -DCMAKE_BUILD_TYPE=Release -S . -B build  # configure

For SCIP, add ``-DENABLE_SCIP`` to the configure step:

    cmake -DENABLE_SCIP -DCMAKE_BUILD_TYPE=Release -S . -B build # configure

If SCIP solver has been installed to the directory other than ``/usr`` or ``/usr/local``, specify the option ``-DCMAKE_MODULE_PATH=/path/to/scip/lib/cmake`` in the configure step.

For HiGHS, add ``-DENABLE_HIGHS`` to the configure step:

    cmake -DENABLE_HIGHS -DCMAKE_BUILD_TYPE=Release -S . -B build # configure

If HiGHS solver has been installed to the directory other than ``/usr`` or ``/usr/local``, specify the option ``-DHiGHS_ROOT=/path/to/highs`` in the configure step.

Usage
-----

RactIP can take two FASTA formatted RNA sequences as input, then
predict their joint secondary structures.

	Usage: ractip [OPTIONS]... [FILES]...

	-h, --help                 Print help and exit
	    --full-help            Print help, including hidden options, and exit
	-V, --version              Print version and exit
	-a, --alpha=FLOAT          weight for hybridization  (default=`0.7')
	-b, --beta=FLOAT           weight for accessibility  (default=`0.0')
  	-t, --fold-th=FLOAT        Threshold for base-pairing probabilities
                               (default=`0.5')
	-u, --hybridize-th=FLOAT   Threshold for hybridazation probabilities
                               (default=`0.1')
	-s, --acc-th=FLOAT         Threshold for accessible probabilities
                               (default=`0.003')
	    --acc-max              optimize for accessibility instead of internal
                               secondary structures  (default=off)
	    --acc-max-ss           additional prediction of interanal secondary
                               structures  (default=off)
	    --acc-num=INT          the number of accessible regions (0=unlimited)
                               (default=`1')
	    --max-w=INT            Maximum length of accessible regions
                               (default=`15')
	    --min-w=INT            Minimum length of accessible regions
                               (default=`5')
	    --zscore=INT           Calculate z-score via dishuffling (0=no shuffling,
                               1=1st seq only, 2=2nd seq only, or 12=both)
                               (default=`0')
	    --num-shuffling=INT    The number of shuffling  (default=`1000')
	    --seed=INT             Seed for random number generator  (default=`0')
	-c, --use-constraint       Use structure constraints  (default=off)
	    --force-constraint     Enforce structure constraints  (default=off)
	    --allow-isolated       Allow isolated base-pairs  (default=off)
	-e, --show-energy          calculate the free energy of the predicted joint
                               structure  (default=off)
	-P, --param-file=FILENAME  Read the energy parameter file for Vienna RNA
                               package

	% ractip DIS.fa DIS.fa
	>DIS
	CUCGGCUUGCUGAGGUGCACACAGCAAGAGGCGAG
	((((.(((((((..[[[[[[.)))))))...))))
	>DIS
	CUCGGCUUGCUGAGGUGCACACAGCAAGAGGCGAG
	((((.(((((((..]]]]]].)))))))...))))

The parentheses '( )' and the brackets '[ ]' indicate the predicted
internal base-pairs and external base-pairs (interactions),
respectively. 

### Run with Docker

    docker build . -t ractip
    docker run -it --rm -v $(pwd):$(pwd) -w $(pwd) ractip ractip DIS.fa DIS.fa

### Run with the web server

The RactIP web server is available [HERE](http://ws.sato-lab.org/rtips/ractip/).

References
----------

* Kato, Y., Sato, K., Hamada, M., Watanabe, Y., Asai, K., Akutsu, T.:
  RactIP: fast and accurate prediction of RNA-RNA interaction using
  integer programming. *Bioinformatics*, 26(18):i460-i466, 2010. [[Link]](https://academic.oup.com/bioinformatics/article/26/18/i460/205351)
* Kato, Y., Mori, T., Sato, K., Maegawa, S., Hosokawa, H., Akutsu, T.:
  An accessibility-incorporated method for accurate prediction of RNA–RNA 
  interactions from sequence data, *Bioinformatics*, 33(2):202-209, 2017. [[Link]](https://academic.oup.com/bioinformatics/article/33/2/202/2525711)
