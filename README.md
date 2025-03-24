# About

`Slope-SpaM` estimates phylogenetic distances between two DNA sequences. More precisely, for a set of input sequences in FASTA format, it estimates all pairwise Jukes-Cantor distances. That is, for each pair of sequences it estimates the average number of substitutions per sequence positions that have occurred since the sequences have evolved from their last common ancestor. The output of the program is a matrix with these pairwise distances in phylip format.

To estimate the phylogenetic distance between two DNA sequences, Slope-SpaM calculates the number of word matches for two different word lengths k_min and k_max. Instead of words, so-called spaced-words may be used, i.e. words containing wildcard characters at certain positions specified by the user. To this end, a binary pattern representing "match positions" ('1') and "don't-care positions" ('0') has to be provided as input to the program.

If spaced words are used, two binary patterns with k_min match positions and k_max match positions are used by shortening the input pattern specified by the user (i.e. by using suitable prefixes of this pattern). 

For more details, check our [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0228070).

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
-t, --threads | Maximum number of threads used.
-a, --as-reads | Use sequences from multi-fasta files as reads of one sequence instead of multiple sequences.
-h, --help | Show help.
input files... | Fasta files or directories of fasta files to process. If more than one input file is given, the option output is ignored. Instead for each .fasta file a .dmat file and for each directory a .dir.dmat file with the same name is created. If --as-reads is not specified, one output matrix is created for each input file. If --as-reads is specified, one output matrix using each input file as a single sequence is created.

## License

Copyright © 2019 - Alexander Linne

License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

## Contact

In case of bugs or unexpected errors don't hesitate to send me a mail: sophie.roehling@stud.uni-goettingen.de
