= RactIP for predicting RNA-RNA interaction using integer programming

== Requirements

* Gurobi Optimizer (>=2.0)
  ((<URL:http://www.gurobi.com/))
* Boost C++ Library (>=1.42.0)
  ((<URL:http://www.boost.org/>))
* Vienna RNA package (>= 1.8)
  ((<URL:http://www.tbi.univie.ac.at/~ivo/RNA/>))

== Install

 ./configure --with-vienna-rna=/path/to/vienna-rna --with-gurobi=/path/to/gurobi
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
