# About

`Slope-SpaM` calculates the number of spaced word matches between pairs of sequences. The matches have *k* identical characters as indicated by a binary pattern where '1' denotes a match position. The number of spaced word matches is observed for a range of *k*-values [k<sub>min</sub> ... k<sub>max</sub>] where k<sub>max</sub> is the weight of the pattern passed to `Slope-SpaM`. The distance is calculated from the slope of a function based on these observed numbers of spaced word matches. For more details, check our [paper](https://www.biorxiv.org/content/10.1101/527515v1).
  
# Installation

The installation requires `cmake` (>= 3.8.2) and a C++17 compliant C++ compiler. Use the following commands to install `Slope-SpaM`:

	git clone --recursive https://github.com/burkhard-morgenstern/Slope-SpaM
	cd Slope-SpaM
	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make

If the default compiler available does not support C++17, the compiler can be set to e.g. GCC 8 using:

	cmake -DCMAKE_C_COMPILER=gcc-8, -DCMAKE_CXX_COMPILER=g++-8 -DCMAKE_BUILD_TYPE=Release ..

# Example Usages

	./slope-spam -o out.dmat -p 10101001010101111111101010101001111111111000111 multi_fasta_file.fasta
	./slope-spam -o out.dmat -k=11,14 multi_fasta_file.fasta
	./slope-spam --as-reads -o out.dmat reads1.fasta reads2.fasta
	./slope-spam --help
	
# Options

Option | Description
:---: | :---
-o, --output | The output file. Ignored if multiple inputs are given.
-p, --pattern | The binary word pattern used to create wordlists from sequences. May only include '0' and '1' characters.
-k, --kmer-lengths | Comma-separated list of wordlengths to consider when calculating the distance between two sequences.
-a, --as-reads | Comma-separated list of wordlengths to consider when calculating the distance between two sequences.
-h, --help | Show help.
input files... | Fasta files or directories of fasta files to process. If more than one input file is given, the option output is ignored. Instead for each .fasta file a .dmat file and for each directory a .dir.dmat file with the same name is created. If --as-reads is not specified, one output matrix is created for each input file. If --as-reads is specified, one output matrix using each input file as a single sequence is created.

## License

Copyright © 2019 - Alexander Linne

License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

## Contact

In case of bugs or unexpected errors don't hesitate to send me a mail: sophie.roehling@stud.uni-goettingen.de
  
