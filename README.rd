= RactIP for predicting RNA-RNA interaction using integer programming

== Requirements

* Boost C++ Library (>=1.42.0) ((<URL:http://www.boost.org/>))
* Vienna RNA package (>= 1.8) ((<URL:http://www.tbi.univie.ac.at/~ivo/RNA/>))
* GNU Linear Programming Kit (>=4.41) ((<URL:http://www.gnu.org/software/glpk/>))
  or Gurobi Optimizer (>=2.0) ((<URL:http://www.gurobi.com/))
  or ILOG CPLEX (>=12.0) ((<URL:http://http://www-01.ibm.com/software/integration/optimization/cplex/>))

== Install

For GLPK,
 ./configure --with-vienna-rna=/path/to/vienna-rna --with-glpk

For Gurobi, 
 ./configure --with-vienna-rna=/path/to/vienna-rna --with-gurobi

For CPLEX,
 ./configure --with-vienna-rna=/path/to/vienna-rna --with-cplex

You may have to specify the include path and the library path by CPPFLAGS and LDFLAGS like
 env CPPFLAGS='-I/path/to/gurobi/include' LDFLAGS='-L/path/to/gurobi/lib' \
 ./configure --with-vienna-rna=/path/to/vienna-rna --with-gurobi

Then,
 make
 make install

== Usage

(({ractip})) can take two FASTA formatted RNA sequences as input, the
predict their joint secondary structures.

 % ractip: [options] fasta1 fasta2
  -h:        show this message
  -p:        do not use the constraints for interenal pseudoknots
  -a alpha:  weight for hybridation probabilities (default: 0.5)
  -t th_bp:  threshold of base-pairing probabilities (default: 0.5)
  -u th_hy:  threshold of hybridazation probabilities (default: 0.2)
  -m:        use McCaskill model (default: CONTRAfold model)
  -i:        allow isolated base-pairs
  -n n_th:   specify the number of threads (default: 1)

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


== References

* Kato, Y., Sato, K., Hamada, M., Watanabe, Y., Asai, K., Akutsu, T.:
  RactIP: fast and accurate prediction of RNA-RNA interaction using
  integer programming. submitted.
