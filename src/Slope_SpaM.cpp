/* Slope Spam ist eine Anwendung um die paarweisen Jukes Cantor Distanzen mehrerer (DNA)-Sequenzen zuberechnen und in einer Phylip-Distanzmatrix auszugeben.
 * Es ist eine allignemtfreie Anwendung, welche auf der Anzahl von (spaced) word matchen beruht.*/

// Importe
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <stdio.h>
#include <map>
#include <math.h>
#include "args.hpp"


// berechnet kmax, das pattern weight
int char_count (std::string pattern){
	int i = 0;
	int kmax = 0;
	int lenPat = pattern.length();
	for  (int i = 0; i < lenPat; i++){
		if (pattern[i] == '1') {
			kmax++;
			}
		}
	return kmax;
}

void create_spaced_words(std::vector<std::vector<int>> & dest, std::string & pattern, std::vector<std::string> & Seqs, size_t seq_num, size_t length)
{
	for (int iter=1; iter<length-pattern.size()+2; iter++){
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
				else if (substring[k] == ' '){
				} 
				else{
					std::cout << "Your sequences also contain other characters than 'A', 'C', 'G' and 'T'. Please checḱ your Input Data."<< '\n';
					std::cout << substring[k] << '\n';
					exit (EXIT_FAILURE);
				}					
			}
		}
		vecword.emplace_back(seq_num);
		dest.emplace_back(vecword);
	}
}

std::pair<float, float> calculate_distance(int (*matche)[2], size_t length1, size_t length2, size_t k_range)
{
	double q = 0.25;
	std::vector <std::vector<long double>> werte;
	for (int w=0; w<k_range; w++){
		if (matche[w][0] != 0){
			long double x = pow(q, matche[w][0]) * length1 * length2;
			long double y = log (matche[w][1] -x);
			std::vector<long double> wertepaar;
			if (isnan(y) == false){
				wertepaar.emplace_back(matche[w][0]);
				wertepaar.emplace_back(y);
			}
			werte.emplace_back(wertepaar);
		}
	}
	double m = 0; 
	double x_mean =0;
	double y_mean=0;
	double x_sum=0;
	double y_sum=0;
	for (int i=0; i<werte.size(); i++){
		x_sum = x_sum + werte[i][0];
	}
	for (int i=0; i<werte.size(); i++){
		y_sum = y_sum + werte[i][1];
	}	
	x_mean = double(1.0/werte.size()) * x_sum;
	y_mean = double(1.0/werte.size()) *y_sum;
	double Zaehler = 0;
	double Nenner = 0;
	for (int i=0; i<werte.size(); i++){
		Zaehler = Zaehler + (werte[i][0]-x_mean)*(werte[i][1]-y_mean);
		Nenner = Nenner + (werte[i][0]-x_mean)*(werte[i][0]-x_mean);
	}
	m = Zaehler/Nenner;
	double p = exp(m);
	double d = -(3.0/4.0) * log(1 - (4.0/3.0) * (1-p));
	return {p, d};
}

std::vector<std::vector<double>> calculate_distance_matrix(std::vector<std::string> & Seqs, std::string & pattern, std::vector<size_t> & Lens, size_t kmax)
{
	std::vector<std::vector<double>> distance(Seqs.size(), std::vector<double>(Seqs.size(), 0));
	// double distance [Seqs.size()] [Seqs.size()] = {0.00000};
	for (int i=0; i<Seqs.size()-1; i++){
		std::cout<<"sequence number : "<<i<<'\n';
		std::vector <std::vector<int>> wordlist1;
		create_spaced_words(wordlist1, pattern, Seqs, i, Lens[i]);
		for (int j=i+1; j<Seqs.size(); j++){
			std::vector <std::vector<int>> wordlistges(wordlist1.begin(), wordlist1.end());
			create_spaced_words(wordlistges, pattern, Seqs, j, Lens[j]);
			std::sort(wordlistges.begin(), wordlistges.end()); 
			int Lmax;
			if (Lens[i]>Lens[j]){
				Lmax = Lens[i];
			}
			else {
				Lmax=Lens[j];
			}
			int kmin =  log (2 * Lmax ) / 0.87 + 1;
			std::map <int, int> matchZahl;
			int matche [kmax-kmin][2]={0};
			for (int k=kmin; k<kmax; k++){
				for (int z=1; z<wordlistges.size(); z++){
					if (z==1){ //erstes paar
						if (std::equal(wordlistges[z-1].begin(), wordlistges[z-1].begin()+k+1, wordlistges[z].begin())){
							matchZahl[wordlistges[z-1][wordlistges[z-1].size()-1]]++;
							matchZahl[wordlistges[z][wordlistges[z].size()-1]]++;
						}
					}
					else { 
						if (std::equal(wordlistges[z-1].begin(), wordlistges[z-1].begin()+k+1, wordlistges[z].begin())){
							if (matchZahl.size()>0){//Match mit vorherigen
								matchZahl[wordlistges[z][wordlistges[z].size()-1]]++;							
							}
							else{
								matchZahl[wordlistges[z-1][wordlistges[z-1].size()-1]]++;
								matchZahl[wordlistges[z][wordlistges[z].size()-1]]++;
							}
						}
						else{//Mismatch
							if (matchZahl.size()>1){
								std::vector<int> v;
								for(std::map<int, int>::iterator iter = matchZahl.begin(); iter != matchZahl.end(); ++iter){
									v.push_back( iter->second );
								}
								int match = v[0]*v[1];
								matche [k-kmin][0] = k+1;
								matche [k-kmin][1] = matche [k-kmin][1] + match;
								matchZahl.clear();	
							}
							else{
								matchZahl.clear();
							}
						}
					}
					if (z==wordlistges.size()-1){//ggf leztes Match eintragen
						if (matchZahl.size()>1){
							std::vector<int> v;
							for(std::map<int, int>::iterator iter = matchZahl.begin(); iter != matchZahl.end(); ++iter){
								v.push_back( iter->second );
							}
							int match = v[0]*v[1];
							matche [k-kmin][0] = k+1;
							matche [k-kmin][1] = matche [k-kmin][1] + match;
							matchZahl.clear();	
						}
						else{
							matchZahl.clear();
						}
					}
				}
			} 
			//y werte berrechnen
			auto values = calculate_distance(matche, Lens[i], Lens[j], kmax - kmin);
			std::cout  << "match probability p : " << values.first << " Jukes-Cantor distance d : " << values.second <<'\n' ;
			distance [i][j] = values.second;
			distance [j][i] = values.second; 
		}
	}
	return distance;
}

int main (int argc, char** argv){
	args::ArgumentParser parser("Slope-SpaM");
	args::ValueFlag<std::string> infile(parser, "Input file", "Specify input file name", {'i', "input"}, args::Options::Required);
	args::ValueFlag<std::string> outfile(parser, "Output file", "Specify output file name", {'o', "output"}, "out.dmat");
	args::ValueFlag<std::string> pattern_flag(parser, "Pattern", "Use this pattern if specified", {'p', "pattern"}, "111111111111111111111111111111111111");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        return 0;
    }
    catch (args::Error e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

	//Einlesen der Sequenzen aus der Fasta-Datei in zwei Vektoren
	std::ifstream file (infile.Get());
	if(file.is_open() == false)
	{
		std::cerr << "File could not be opened: " << infile.Get() << std::endl;
		return -1;
	}
	std::vector <std::string> SeqKeys; // Vektor für die Namen der Sequenzen
	std::vector <std::string> Seqs; // Vektor für die Sequenzen als ein String
	std::string SeqK;
	std::string SeqKey;
	std::string Seq;
	std::string line;
	while (getline(file, line)){
		if (line [0] == '>') {
			if (!Seq.empty()){
				Seqs.emplace_back(Seq);
			}
			SeqK = line;
			Seq = '\0';
			SeqKey = SeqK.substr(1);
			SeqKeys.emplace_back(SeqKey);
		}
		else {
			Seq = Seq + line;
		}
	}
	
	Seqs.emplace_back(Seq);
	std::cout << "Data imported" << '\n'; // Die Sequenzen und ihre Namen wurden in den entsprechenden Vektoren gespeichert
	// Berechnen der Längen aller Sequenzen und Speichern in einem Vektor
	std::vector <size_t> Lens; // Vektor für die Längen der Sequenzen
	size_t len;
	for (int i=0; i<Seqs.size(); i++){
		len = Seqs[i].length()-1;
		Lens.emplace_back(len);
	}
	// Verarbeiten des Patterns
	std::string pattern = pattern_flag.Get();
	int lenPat = pattern.length();  // Länge des Pattern
	int kmax; // pattern weight, Anzahl der Matchpositionen im Pattern
	kmax = char_count(pattern);
	int Lm  = *std::max_element(Lens.begin(), Lens.end());
	int km =  log (2 * Lm ) / 0.87 + 1;
	if (kmax-km < 10){
		std::cout << "Your pattern weight is too small for your sequences, there should be a least 10 value calculable. Please use a pattern with a higher pattern weight."<< '\n';
		exit (EXIT_FAILURE);
	}
	int SeqZahl = Lens.size();	// Anzahl der Sequenzen	
	// Erstellen der Outputdatei und eintragen der Anzahl der Sequenzen
	std::ofstream f (outfile.Get());
	f << SeqZahl << '\n';
	f.close();
	
	auto distance = calculate_distance_matrix(Seqs, pattern, Lens, kmax);
	std::ofstream ff (outfile.Get(), std::ios::app);
	for (int i=0; i<SeqZahl; i++){
		ff << SeqKeys [i] << '\t';
		for (int k=0; k<SeqZahl; k++){
			ff << distance[i][k] <<'\t';
		}
		ff <<'\n';
	}
	ff.close();
}