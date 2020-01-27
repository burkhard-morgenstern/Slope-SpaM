# Slope-SpaM  &nbsp;
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/ebcd8da0747b48da84693965306a4ac8)](https://www.codacy.com/manual/smortezah/Slope-SpaM?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=burkhard-morgenstern/Slope-SpaM&amp;utm_campaign=Badge_Grade)
![GitHub repo size](https://img.shields.io/github/repo-size/burkhard-morgenstern/Slope-SpaM)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)

`Slope-SpaM` estimates phylogenetic distances between two DNA sequences. More precisely, for a set of input sequuences in FASTA format, it estimates all pairwise Jukes-Cantor distances. That is, for each pair of sequences it estimates the average number of substitutions per sequence positions that have occurred since the sequences have evolved from their last common anchestor. The output of the program is a matrix with these pairwise distances in phylip format.

To estimate the phylogenetic distance between two DNA sequences, Slope-SpaM calculates the number of word matches for two different word lengths k_min and k_max. Instead of words, so-called spaced-words may be used, i.e. words containing wildcard characters at certain positions specified by the user. To this end, a binary pattern representing "match positions" ('1') and "don't-care positions" ('0') has to be provided as input to the program.

If spaced words are used, two binary patterns with k_min match positions and k_max match positions are used by shortening the input pattern specified by the user (i.e. by using suitable prefixes of this pattern). 

For more details, check our [paper](https://www.biorxiv.org/content/10.1101/527515v1).

## Installation
The installation requires CMake >= 3.8.2 and a C++17 compliant compiler.

### Linux / macOS
CMake v3.15.5 can be downloaded [here](https://github.com/Kitware/CMake/releases/download/v3.15.5/cmake-3.15.5-Linux-x86_64.sh) for Linux, and [here](https://github.com/Kitware/CMake/releases/download/v3.15.5/cmake-3.15.5-Darwin-x86_64.dmg) for macOS. Use the following commands to install `Slope-SpaM`; the binaries will be provided in the `bin/` folder:

	git clone --recursive https://github.com/burkhard-morgenstern/Slope-SpaM
	cd Slope-SpaM
	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make
	mkdir ../bin
	mv slope-spam slope-spam-test generator simulate ../bin

If the default compiler available does not support C++17, the compiler can be set to e.g. GCC 8 using:

	cmake -DCMAKE_C_COMPILER=gcc-8, -DCMAKE_CXX_COMPILER=g++-8 -DCMAKE_BUILD_TYPE=Release ..

### Windows 10
After installing `wsl` (Windows Subsystem for Linux), you can use the same commands as for Linux/macOS.

## Example Usages

	./slope-spam -o out.dmat -p 10101001010101111111101010101001111111111000111 multi_fasta_file.fasta
	./slope-spam -o out.dmat -k=11,14 multi_fasta_file.fasta
	./slope-spam --as-reads -o out.dmat reads1.fasta reads2.fasta
	./slope-spam --help

## Options
| Option                        | Description |
| ---                           | ---         |
| -o, --output <img width=600/> | The output file. Ignored if multiple inputs are given. <img width=600/> |
| -p, --pattern                 | The binary word pattern used to create wordlists from sequences. May only include '0' and '1' characters. |
| -k, --kmer-lengths            | Comma-separated list of wordlengths to consider when calculating the distance between two sequences. |
| -a, --as-reads                | Comma-separated list of wordlengths to consider when calculating the distance between two sequences. |
| -h, --help                    | Show help. |
| input files...                | FASTA files or directories of FASTA files to process. If more than one input file is given, the option output is ignored. Instead for each .fasta file a .dmat file and for each directory a .dir.dmat file with the same name is created. If --as-reads is not specified, one output matrix is created for each input file. If --as-reads is specified, one output matrix using each input file as a single sequence is created. |

### License
Copyright © 2019 Alexander Linne

License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

### Contact
In case of bugs or unexpected errors don't hesitate to send me a mail: sophie.roehling@stud.uni-goettingen.de
