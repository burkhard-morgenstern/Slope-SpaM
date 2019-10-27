# About

  

`Slope-SpaM` estimates phylogenetic distances between two DNA sequences. More precisely, for a set of input sequuences in FASTA format, it estimates all pairwise Jukes-Cantor distances. That is, for each pair of sequences it estimates the average number of substitutions per sequence positions that have occurred since the sequences have evolved from their last common anchestor. The output of the program is a matrix with these pairwise distances in phylip format.

To estimate the phylogenetic distance between two DNA sequences, Slope-SpaM calculates the number of word matches for two different word lengths k_min and k_max. Instead of words, so-called spaced-words may be used, i.e. words containing wildcard characters at certain positions specified by the user. To this end, a binary pattern representing "match positions" ('1') and "don't-care positions" ('0') has to be provided as input to the program.

If spaced words are used, two binary patterns with k_min match positions and k_min match positions are used by shortening the input pattern specified by the user (i.e. by using suitable prefixes of this pattern). 

For more details, check our [paper](https://www.biorxiv.org/content/10.1101/527515v1).
  

# Installation and Usage

  
The installation  requires `cmake`. Use the following commands to install `Slope-SpaM`:

	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make
	cd ..
Example usage:

	./slope-spam -i <path_to_fasta_file> -o distance_matrix -p 10101001010101111111101010101001111111111000111
	./slope-spam -i <path_to_fasta_file> -o distance_matrix -p <path_to_pattern_file>
	
# Options

  

Option | Description
:---: | :---
-i | Input file in FASTA format
-o | Output file in phylip format (distance matrix)
-p | Binary pattern (either directly or a filename containing the pattern in the first line)
  

## License

  

Copyright © 2019 - Sophie Roehling

License GPLv3+: GNU GPL version 3 or later.

  

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

## Contact

  

In case of bugs or unexpected errors don't hesitate to send me a mail: sophie.roehling@stud.uni-goettingen.de
  
