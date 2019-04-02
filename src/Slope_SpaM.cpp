/* Slope Spam ist eine Anwendung um die paarweisen Jukes Cantor Distanzen mehrerer (DNA)-Sequenzen zuberechnen und in einer Phylip-Distanzmatrix auszugeben.
 * Es ist eine allignemtfreie Anwendung, welche auf der Anzahl von (spaced) word matchen beruht.*/

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <map>
#include <cmath>
#include "args.hpp"

void parse_options(args::ArgumentParser & parser, int argc, char ** argv)
{
	try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        exit(0);
    }
    catch (args::Error e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(-1);
    }
}

void read_sequences(std::string & file_name, std::vector<std::string> & sequences, std::vector<std::string> & sequence_names, std::vector<size_t> & lengths)
{
    std::ifstream infile(file_name);
    if ( infile.is_open() == false )
    {
        std::cerr << "File " + file_name + " could not be opened!";
		exit(-1);
    }
    std::string line;
    std::string header;
    std::getline( infile, line, '>' );
    while ( !infile.eof() )
    {
        std::getline( infile, header );
        std::getline( infile, line, '>' );
		line.erase(std::remove_if(line.begin(), line.end(), [](const char c ){return isspace(c);}), line.end());
		sequence_names.emplace_back(header);
        sequences.emplace_back(line);
		lengths.push_back(line.size());
	}
	std::cout << sequences.size() << " sequences have been read." << std::endl;
}

int check_pattern(std::string & pattern, std::vector<size_t> & Lens)
{
	int lenPat = pattern.length();
	int kmax = std::count_if(pattern.begin(), pattern.end(), [](const char c){return c == '1';});
	int Lm  = *std::max_element(Lens.begin(), Lens.end());
	int km =  log (2 * Lm ) / 0.87 + 1;
	if (kmax-km < 10){
		// TODO: was does this even mean?
		std::cout << "Your pattern weight is too small for your sequences, there should be a least 10 value calculable. Please use a pattern with a higher pattern weight."<< '\n';
		exit(-1);
	}
	return kmax;
}

void create_spaced_words(std::vector<std::vector<int>> & dest, std::string & pattern, std::vector<std::string> & Seqs, size_t seq_num, size_t length)
{
	for (int iter=0; iter<length-pattern.size()+2; iter++){
		auto substring = Seqs[seq_num].begin() + iter;
		std::vector <int> vecword;
		for (int k=0; k<pattern.size(); k++){
			if (pattern[k] == '1'){
				if (substring[k] == 'A'){
					vecword.push_back(1);
				}
				else if (substring[k] == 'C'){
					vecword.push_back(2);
				}
				else if (substring[k] == 'G'){
					vecword.push_back(3);
				}
				else if (substring[k] == 'T'){
					vecword.push_back(4);
				}
				else{
					break;
				}					
			}
		}
		if(vecword.size() < pattern.size())
		{
			continue; // illegal characters
		}
		vecword.emplace_back(seq_num);
		dest.emplace_back(vecword);
	}
}

std::pair<float, float> calculate_distance(std::vector<size_t> & matches, size_t length1, size_t length2, size_t kmin)
{
	double q = 0.25;
	std::vector <std::pair<long double, long double>> werte;
	for (int w=0; w<matches.size(); w++){
		if (matches[w] != 0){
			long double x = pow(q, kmin + w + 1) * length1 * length2; // TODO:again, why + 1
			long double y = log (matches[w] - x);
			if (std::isnan(y) == false){
				werte.emplace_back(kmin + w + 1, y);
			}
		}
	}
	double m = 0; 
	double x_mean =0;
	double y_mean=0;
	double x_sum=0;
	double y_sum=0;
	for (int i=0; i<werte.size(); i++){
		x_sum += werte[i].first;
	}
	for (int i=0; i<werte.size(); i++){
		y_sum += werte[i].second;
	}	
	x_mean = double(1.0/werte.size()) * x_sum;
	y_mean = double(1.0/werte.size()) *y_sum;
	double Zaehler = 0;
	double Nenner = 0;
	for (int i=0; i<werte.size(); i++){
		Zaehler += (werte[i].first-x_mean)*(werte[i].second-y_mean);
		Nenner += (werte[i].first-x_mean)*(werte[i].first-x_mean);
	}
	m = Zaehler/Nenner;
	double p = exp(m);
	double d = -(3.0/4.0) * log(1 - (4.0/3.0) * (1-p));
	return {p, d};
}

std::vector<std::vector<double>> calculate_distance_matrix(std::vector<std::string> & Seqs, std::string & pattern, std::vector<size_t> & Lens, size_t kmax)
{
	std::vector<std::vector<double>> distance(Seqs.size(), std::vector<double>(Seqs.size(), 0));
	for (int i=0; i<Seqs.size()-1; i++)
	{
		std::cout<<"sequence number : "<<i<<'\n';
		std::vector <std::vector<int>> wordlist1;
		create_spaced_words(wordlist1, pattern, Seqs, i, Lens[i]);
		for (int j=i+1; j<Seqs.size(); j++)
		{
			std::vector <std::vector<int>> wordlistges(wordlist1.begin(), wordlist1.end());
			create_spaced_words(wordlistges, pattern, Seqs, j, Lens[j]);
			std::sort(wordlistges.begin(), wordlistges.end()); 
			int Lmax = std::max(Lens[i], Lens[j]);
			int kmin =  log (2 * Lmax ) / 0.87 + 1;
			std::vector<size_t> matches(kmax - kmin);
			for (int k = kmin; k < kmax; k++) // TODO: kmax not used?
			{
				for(auto it = wordlistges.begin(); it != wordlistges.end();)
				{
					auto next = std::find_if_not(it, wordlistges.end(), [=](const std::vector<int> & word)
																		{return std::equal(word.begin(), word.begin() + k + 1, it->begin());});
					int seq_i_count = std::count_if(it, next, [=](std::vector<int> & word){return word.back() == i;});
					int seq_j_count = std::distance(it, next) - seq_i_count;
					// matche[k - kmin][0] = k + 1; // TODO: why was there k+1 here?
					matches[k - kmin] += seq_i_count * seq_j_count;
					it = next;
				}
			} 
			//y werte berrechnen
			auto values = calculate_distance(matches, Lens[i], Lens[j], kmin);
			std::cout  << "match probability p : " << values.first << " Jukes-Cantor distance d : " << values.second <<'\n' ;
			distance [i][j] = values.second;
			distance [j][i] = values.second; 
		}
	}
	return distance;
}

void print_distance_matrix(std::string outfile, size_t num_sequences, std::vector<std::string> & SeqKeys, std::vector<std::vector<double>> & distance)
{
	std::ofstream ff (outfile, std::ios::app);
	ff << num_sequences << std::endl;
	for (int i=0; i<num_sequences; i++){
		ff << SeqKeys [i] << '\t';
		for (int k=0; k<num_sequences; k++){
			ff << distance[i][k] <<'\t';
		}
		ff << std::endl;
	}
	ff.close();
}

int main (int argc, char** argv){
	args::ArgumentParser parser("Slope-SpaM");
	args::ValueFlag<std::string> infile(parser, "Input file", "Specify input file name", {'i', "input"}, args::Options::Required);
	args::ValueFlag<std::string> outfile(parser, "Output file", "Specify output file name", {'o', "output"}, "out.dmat");
	args::ValueFlag<std::string> pattern_flag(parser, "Pattern", "Use this pattern if specified", {'p', "pattern"}, "111111111111111111111111111111111111");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	parse_options(parser, argc, argv);

	std::vector <std::string> SeqKeys;
	std::vector <std::string> Seqs;
	std::vector <size_t> Lens;
	read_sequences(infile.Get(), Seqs, SeqKeys, Lens);

	std::string pattern = pattern_flag.Get();
	int kmax = check_pattern(pattern, Lens);

	auto distance = calculate_distance_matrix(Seqs, pattern, Lens, kmax);
	print_distance_matrix(outfile.Get(), Seqs.size(), SeqKeys, distance);
}