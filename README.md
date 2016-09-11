# Alchemy 2 - Inference and Learning in Markov Logic
This is a fork of the code located at https://code.google.com/p/alchemy-2/.


## Description

Alchemy is a software package providing a series of algorithms for statistical relational learning and probabilistic logic inference, based on the Markov logic representation. Alchemy allows you to easily develop a wide range of AI applications, including:

* Collective classification
* Link prediction
* Entity resolution
* Social network modeling
* Information extraction

More info at http://alchemy.cs.washington.edu/


## Code

* src/ contains source code and a makefile.
* doc/ contains a change log, and a manual in PDF, PostScript and html formats.
* exdata/ contains a simple example of Alchemy input files.
* bin/ is used to contain compiled executables.


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

## Building
```
$ cd path/to/alchemy-2/
$ cd src
$ make depend
$ make
```

Note: This fork of http://code.google.com/p/alchemy-2 has been updated to compile properly on a Mac following instructions from http://alchemy.cs.washington.edu/requirements.html


## License

By using Alchemy, you agree to accept the license agreement in LICENSE.md
