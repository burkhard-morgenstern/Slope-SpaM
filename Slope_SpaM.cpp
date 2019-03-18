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


// Das Eigentliche Programm
int main (int argc, char** argv){
	// Abrage alle nötigen Eingaben
	// Eingabe der Fasta-Datei der Sequenzen
	std::string filename;
	std::string yn;
	do {
		std::cout << "Please enter the name of the Fasta file containing the sequences you want to compare." << '\n';
		std::cin >> filename;
		std::cout << "You entered " << filename << " as the filename. Is that correct? (y/n)" << '\n';
		std::cin >> yn;
	} while( yn == "n" );
	// Eingabe des Patterns
	std::string pattern;
	do {
		std::cout << "Please enter the pattern, consisting of ones, for match positions, and zeroes, for do not care positions." << '\n';
		std::cin >> pattern;
		std::cout << "You entered " << pattern << " as pattern. Is that correct? (y/n)" << '\n';
		std::cin >> yn;
	} while( yn == "n" );
	// Eingabe des Names für die Output Datei
	std::string outfile;
	do {
		std::cout << "Please enter a name for your output file, which will contain the distance matrix."<< '\n';
	std::cin >> outfile;
	std::cout << "You entered " << outfile << " as name for the output file. Is that correct? (y/n)" << '\n';
		std::cin >> yn;
	} while( yn == "n" );	
	//Einlesen der Sequenzen aus der Fasta-Datei in zwei Vektoren
	std::ifstream file (filename);
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
	std::vector <int> Lens; // Vektor für die Längen der Sequenzen
	int len;
	for (int i=0; i<Seqs.size(); i++){
		len = Seqs[i].length()-1;
		Lens.emplace_back(len);
	}
	// Verarbeiten des Patterns
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
	std::ofstream f (outfile);
	f << SeqZahl << '\n';
	f.close();
	double distance [SeqZahl] [SeqZahl] = {0.00000};
	for (int i=0; i<Seqs.size()-1; i++){
		std::cout<<"sequence number : "<<i<<'\n';
		std::vector <std::vector<int>> wordlist1;
		for (int iter=1; iter<Lens[i]-lenPat+2; iter++){
			std::string substring = Seqs[i].substr(iter, lenPat);
			std::string word;
			std::vector <int> vecword;
			for (int k=0; k<lenPat; k++){
				if (pattern[k] == '1'){
					if (substring[k] == 'A'){
						word.append("1");
					}
					else if (substring[k] == 'C'){
						word.append("2");
					}
					else if (substring[k] == 'G'){
						word.append("3");
					}
					else if (substring[k] == 'T'){
						word.append("4");
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
			for (int l=0; l<word.length(); l+=1){
				std::string tw = word.substr(l, 1);
				int zahlword =  std::atoi (tw.c_str());
				vecword.emplace_back(zahlword);
			}
			vecword.emplace_back(i);
			wordlist1.emplace_back(vecword);
		}
		for (int j=i+1; j<Seqs.size(); j++){
			std::vector <std::vector<int>> wordlist2;
			for (int iter=1; iter<Lens[j]-lenPat+2; iter++){
				std::string substring = Seqs[j].substr(iter, lenPat);
				std::string word;
				std::vector <int> vecword;
				for (int k=0; k<lenPat; k++){
					if (pattern[k] == '1'){
						if (substring[k] == 'A'){
							word.append("1");
						}
						else if (substring[k] == 'C'){
							word.append("2");
						}
						else if (substring[k] == 'G'){
							word.append("3");
						}
						else if (substring[k] == 'T'){
							word.append("4");
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
				for (int l=0; l<word.length(); l+=1){
					std::string tw = word.substr(l, 1);
					int zahlword =  std::atoi (tw.c_str());
					vecword.emplace_back(zahlword);
				}
				vecword.emplace_back(j);
				wordlist2.emplace_back(vecword);
			}
			std::vector <std::vector<int>> wordlistges;
			wordlistges.insert (wordlistges.end(), wordlist1.begin(), wordlist1.end());
			wordlistges.insert (wordlistges.end(), wordlist2.begin(), wordlist2.end());
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
							matchZahl[wordlistges[z-1][wordlistges[z-1].size()-1]] = matchZahl[wordlistges[z-1][wordlistges[z-1].size()-1]]+1;
							matchZahl[wordlistges[z][wordlistges[z].size()-1]] = matchZahl[wordlistges[z][wordlistges[z].size()-1]]+1;
						}
					}
					else { 
						if (std::equal(wordlistges[z-1].begin(), wordlistges[z-1].begin()+k+1, wordlistges[z].begin())){
							if (matchZahl.size()>0){//Match mit vorherigen
								matchZahl[wordlistges[z][wordlistges[z].size()-1]] = matchZahl[wordlistges[z][wordlistges[z].size()-1]]+1;								
							}
							else{
								matchZahl[wordlistges[z-1][wordlistges[z-1].size()-1]] = matchZahl[wordlistges[z-1][wordlistges[z-1].size()-1]]+1;
								matchZahl[wordlistges[z][wordlistges[z].size()-1]] = matchZahl[wordlistges[z][wordlistges[z].size()-1]]+1;
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
		double q = 0.25;
		std::vector <std::vector<long double>> werte;
		for (int w=0; w<kmax-kmin; w++){
			if (matche[w][0] != 0){
			long double x = pow(q, matche[w][0]) * Lens[i] * Lens[j];
			long double y = log (matche[w][1] -x);
			std::vector<long double> wertepaar;
			if (isnan(y)){
			}
			else{
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
		double Zaehler;
		double Nenner;
		for (int i=0; i<werte.size(); i++){
			Zaehler = Zaehler + (werte[i][0]-x_mean)*(werte[i][1]-y_mean);
			Nenner = Nenner + (werte[i][0]-x_mean)*(werte[i][0]-x_mean);
		}
		m = Zaehler/Nenner;
		double p = exp(m);
		double d = -(3.0/4.0) * log(1 - (4.0/3.0) * (1-p));
		std::cout  << "match probability p : " << p << "Jukes-Cantor distance d : " << d <<'\n' ;
		distance [i][j] = d;
		distance [j][i] = d; 
		wordlist2.clear();
		wordlistges.clear();
}
}
	std::ofstream ff (outfile, std::ios::app);
	for (int i=0; i<SeqZahl; i++){
		ff << SeqKeys [i] << '\t';
		for (int k=0; k<SeqZahl; k++){
			ff << distance[i][k] <<'\t';
		}
		ff <<'\n';
	}
	ff.close();
	std::cout << "finished"<<'\n';
}
	















