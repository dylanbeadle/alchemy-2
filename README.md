# Alchemy 2 - Inference and Learning in Markov Logic
This is a fork of the code located at https://code.google.com/p/alchemy-2/.


## Description

Alchemy is a software package providing a series of algorithms for statistical relational learning and probabilistic logic inference, based on the Markov logic representation. Alchemy allows you to easily develop a wide range of AI applications, including:

* Collective classification
* Link prediction
* Entity resolution
* Social network modeling
* Information extraction

If you are not already familiar with Markov logic, we recommend that you first read the paper [Unifying Logical and Statistical AI](http://www.cs.washington.edu/homes/pedrod/papers/aaai06c.pdf). If you want to understand how lifted inference algorithms operate, read the [Probabilistic Theorem Proving](http://www.hlt.utdallas.edu/~vgogate/papers/uai11-b.pdf) paper.

Alchemy 2.0 includes the following algorithms:

* Discriminative weight learning (Voted Perceptron, Conjugate Gradient, and Newton's Method)
* Generative weight learning
* Structure learning
* propositional MAP/MPE inference (including memory efficient)
* propositional and lazy Probabilistic inference algorithms: MC-SAT, Gibbs Sampling and Simulated Tempering
* Lifted Belief propagation
* Support for native and linked-in functions
* Block inference and learning over variables with mutually exclusive and exhaustive values
* EM (to handle ground atoms with unknown truth values during learning)
* Specification of indivisible formulas (i.e. formulas that should not be broken up into separate clauses)
* Support of continuous features and domains
* Online inference
* Decision Theory
* Probabilistic theorem proving (lifted weighted model counting)
* Lifted importance sampling
* Lifted Gibbs sampling

More info at http://alchemy.cs.washington.edu/

## Code

* ```src/``` contains source code and a makefile.
* ```doc/``` contains a change log, and a manual in PDF, PostScript and html formats.
* ```exdata/``` contains a simple example of Alchemy input files.
* ```bin/``` is used to contain compiled executables.

## Dependencies
* g++ 4.1.2
* Bison 2.3
* Flex 2.5.4
* Perl 5.8.8

You can install perl and gcc using Homebrew on Mac. Bison and Flex must be present already.
```
$ brew tap homebrew/versions
$ brew install gcc49 perl518
```
## Build
Either git-clone or extract the downloaded archive in $PROJECT_HOME
  * ```cd $PROJECT_HOME/src```
  * ```make depend```
  * ```make```


Note: This fork of http://code.google.com/p/alchemy-2 has been updated to compile properly on a Mac following instructions from http://alchemy.cs.washington.edu/requirements.html

## Usage

### Structure learning
Learn the structure of a model given a training database consisting of ground atoms
```
learnstruct -i <input .mln file> -o <output .mln file> -t <training .db file>
```

### Weight learning 
Learn parameters of a model given a training database consisting of ground atoms
```
learnwts -i <input .mln file> -o <output .mln file> -t <training .db file>
```

### *Inference* 
Infer the probability or most likely state of query atoms given a test database consisting of evidence ground atoms
```
infer -i <input .mln file> -r <output file containing inference results> -e <evidence .db file> -q <query atoms (comma-separated with no space)>
```

## Tutorial

Tutorial: https://alchemy.cs.washington.edu/tutorial/tutorial.html
More: http://alchemy.cs.washington.edu/

## License

By using Alchemy, you agree to accept the license agreement in LICENSE.md
