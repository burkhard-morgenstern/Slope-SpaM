# About

  

`Slope-SpaM` calculates the number of spaced word matches between pairs of sequences. The matches have *k* identical characters as indicated by a binary pattern where '1' denotes a match position. The number of spaced word matches is observed for a range of *k*-values [k<sub>min</sub> ... k<sub>max</sub>] where k<sub>max</sub> is the weight of the pattern passed to `Slope-SpaM`. The distance is calculated from the slope of a function based on these observed numbers of spaced word matches. For more details, check our [paper](https://www.biorxiv.org/content/10.1101/527515v1).
  

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
--- | ---
-i | Input file in FASTA format
-o | Output file in phylip format (distance matrix)
-p | Binary pattern (either directly or a filename containing the pattern in the first line)
  

## License

  

Copyright © 2019 - Sophie Roehling

License GPLv3+: GNU GPL version 3 or later.

  

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

## Contact

  

In case of bugs or unexpected errors don't hesitate to send me a mail: sohpie.roehling@stud.uni-goettingen.de)
  
