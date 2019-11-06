/*
 * generator.cpp
 *
 *  Created on: 2 May 2016
 *	Author: salma
 */
//Sequence evolution simulation

#include "generator.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <getopt.h>

using namespace std;

//vector<vector<char>> alignment (2, vector<char>());

void printHelp(){
 std::string help = 
    "\nUsage: ./seqGen [options] "
    "\nOptions:"        
    "\n\t -h: print this help and exit"
    "\n\t -m <double>: substitution rate"
    "\n\t -i <double>: indel rate"
    "\n\t -f <file>: use genome in <file> as root genome, if -f unspecified an artificial sequence is generated as root"
	"\n\t -l <integer>: specify length if no root genome is provided, i.e. if -f <file> is unspecified"
    "\n";
	std::cout << help << std::endl;
}

//gets three arguments: 0.name of the file 1.mutation rate 2.indel rate
//main can be modified to accept more or fewer arguments. we just need to keep in mind argv starts from 1
#ifndef GENERATOR_NO_MAIN
int main(int argc, char *argv[]) {

    clock_t t1,t2;
    //t1=clock();

	//in case of reading ancestral seq from a file this should be commented out:
	//string input(argv[1]);
    string input = "";
	//input argument: mutation rate
	double mRate = 0.1;

	//input argument: indel rate
	double idRate = 0;
	int length = 5000000;
	//input argument: repeat rate
	double rRate = 0; //atof(argv[3]);
	
	//input argument: cutting parameters(in case we want one of descendents with a proportion of ansectral genome)
	//double startToCut = atof(argv[3]);
	//double cutRatio = atof(argv[4]);
	string output_path;
	int option_char;
	 while ((option_char = getopt (argc, argv, "m:i:f:l:r:o:")) != -1){ 
		switch (option_char){  
			case 'm': 
				mRate = atof (optarg); 
				break;
			case 'l': 
				length = atoi (optarg); 
				break;
			case 'i': 
				idRate = atof (optarg); 
				break;
			case 'r': 
				rRate = atof (optarg); 
				break;
			case 'f':
				input = optarg;
				break;
			case 'o':
				output_path = optarg;
				break;
			case 'h': 
				printHelp();
				exit (EXIT_SUCCESS);
				break;
			case '?': 
				printHelp();		
				exit (EXIT_FAILURE);
      	}
	}
	vector<char> sequence;
	//in case of artificial seq 
	if(input.length() == 0)
		sequence = DNA_SeqGenerator (length);
	else
		sequence = readFASTASeqFromFile (input);
	//in case of reading from file
	//vector<char> sequence = readFASTASeqFromFile (input);


	//vector<char> dseq1 = mutate(sequence, mRate/2.0);
	//vector<char> dseq2 = mutate(sequence, mRate/2.0);
	vector<char> dseq1 = sequence;
	vector<char> dseq2 = mutate(sequence, mRate);

	//dseq2 = cutSequence (dseq2, startToCut, cutRatio);
	//cout << "Finished mutations\n";

//	alignment[0] = dseq1;
//	alignment[1] = dseq2;

	//indels can be imposed on one sequence or both(idRate/2)
	if(idRate>0)
		dseq1 = indel (dseq1 , idRate);

	if (rRate > 0) {
		dseq1 = repeat (dseq1, rRate/2);
		dseq2 = repeat (dseq2, rRate/2);
	}

	//cout << "Finished indels\n";

	//writing the descending sequences into a fasta file. name of the file can be changed
	saveFASTA (dseq1, dseq2, output_path);

	//can make an alignment after indels but the parts in other functions that have been commented
	//out need to be added again
//	printAlignment (alignment , "align.txt");

//	double d = Jukes_Cantor(alignment[0], alignment[1]);

//	Print(sequence);
//	Print(dseq1);
//	Print(dseq2);


    //t2=clock();
    //cout << "time = " << (t2-t1)/ CLOCKS_PER_SEC << endl;
//	cout << "Jukes_Cantor distance estimate:\t" << d << endl;

	return 0;
}
#endif

//Generate a random DNA sequence with equal A:T:C:G ratios
vector<char> DNA_SeqGenerator(unsigned long int length)
{
	vector<char> seq (length,'0');

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, 3);

    int freq[4] = {0,0,0,0};
    char bases[4] = {'A' , 'T' , 'C' , 'G'};

    int cur=0;
    for (unsigned long int i=0; i<length; i++)
    {
        cur = dist(mt);
        seq[i] = bases[cur];
        freq[cur]++;
    }

    //cout << "A Random Sequence was Generated: " << "Frequencies are:\t" << freq[0] << "\t" << freq[1] << "\t" << freq[2] << "\t" << freq[3] << endl;

    return seq;
}


//Mutate the ancestral sequence with the specified mutation rate
vector<char> mutate(vector<char> motherSeq, double mutationRate)
{
	vector<char> ds = motherSeq;

	//ratio of ti/tv
	const double ratio = 2.0;

	//random number generation
	random_device rd1;
	mt19937 mt1(rd1());
	uniform_int_distribution<unsigned long int> dist1(0, motherSeq.size()-1);

	//generate random double no to determine substitute base based on ti/tv rate
	random_device rd2;
	mt19937 mt2(rd2());
	uniform_real_distribution<double> dist2(0.0, (2*ratio + 2));


	unsigned long int index = 0;
	double change = 0;

    static const map<char, vector<char>> substitute =
	{
		{'A', {'G', 'T', 'C'}},
		{'T', {'C', 'A', 'G'}},
		{'C', {'T', 'A', 'G'}},
		{'G', {'A', 'T', 'C'}},
		{'N', {'X', 'X', 'X'}}
	};


	for(unsigned long int i=0 ; i < mutationRate* motherSeq.size() ; i++)
	{
		index = dist1(mt1);
		change = dist2(mt2);

		//x = the character to be replaced
		char x = motherSeq [index];

		if (x == 'N')
		{
			i--;
			continue;
		}
		//replace x based on the probabilities indicated by ratio
		if (change < ratio*2)
		{
			ds[index] = substitute.at(x)[0];
		}
		else {if ((change > ratio*2) & (change < (ratio*2 + 1)))
				ds[index] = substitute.at(x)[1];
			else
				ds[index] = substitute.at(x)[2];
		}
	}

	//cout << "\n";

	return ds;
}

//Mutate the ancestral sequence with the specified mutation rate
vector<char> mutate2(vector<char> motherSeq, double mutationRate)
{
	vector<char> ds = motherSeq;

	//ratio of ti/tv
	const double ratio = 2.0;

	//random number generation
	random_device rd1;
	mt19937 mt1(rd1());
	uniform_real_distribution<double> dist1(0, 1);

	//generate random double no to determine substitute base based on ti/tv rate
	random_device rd2;
	mt19937 mt2(rd2());
	uniform_real_distribution<double> dist2(0.0, (2*ratio + 2));

    vector<char> bases = {'A', 'T', 'C', 'G', 'N'};

    vector<vector<char>> substitute (5, vector<char>(3, '0'));
    substitute[0] = {'G', 'T', 'C'};
    substitute[1] = {'C', 'A', 'G'};
    substitute[2] = {'T', 'A', 'G'};
    substitute[3] = {'A', 'T', 'C'};
    substitute[4] = {'X', 'X', 'X'};


	for(unsigned long int i = 0; i < motherSeq.size(); i++)
	{
		auto p = dist1(mt1);
		if (p > mutationRate)
			continue;

		auto change = dist2(mt2);

		//x = the character to be replaced
		auto x = motherSeq[i];

		//pos = position of replaced character in bases vector
		unsigned long int pos = find(bases.begin(), bases.end(), x) - bases.begin();

		if (pos == 4)
		{
			i--;
			continue;
		}
		//replace x based on the probabilities indicated by ratio
		if (change < ratio*2)
		{
			ds[i] = substitute[pos][0];
		}
		else {if ((change > ratio*2) & (change < (ratio*2 + 1)))
				ds[i] = substitute[pos][1];
			else
				ds[i] = substitute[pos][2];
		}
	}

	//cout << "\n";

	return ds;
}

//calculate the distance between two seqs
double Jukes_Cantor(vector<char>& s1, vector<char>& s2){
	int mutationNumber = 0;
	int length = 0;
	for (unsigned int i=0; i<s1.size();i++){
		if ((s1[i]!='-') & (s2[i]!='-'))
			length ++;
		if ( (s1[i]!=s2[i]) & (s1[i]!='-') & (s2[i]!='-') )
			mutationNumber++;
	}
	double p = mutationNumber/((double)length);
	double d = -0.75*(log(1-(4.0/3.0)*p));
	return d;
}

//prints a vector
void Print (vector<char>& v){

  for (unsigned int i=0; i<v.size();i++){
    cout << v[i];
  }
  cout <<endl;
}

//saving the fasta file of the two new sequences
void saveFASTA (vector<char> seq1, vector<char> seq2, string name)
{
	ofstream fastafile;
	fastafile.open (name.c_str());

	fastafile << ">" << "seq1";
	for (unsigned long int i=0; i<seq1.size(); i++){
		if((i%100)==0)
			fastafile << "\n";
		fastafile << seq1[i];
	}

	fastafile <<"\n";

	fastafile << ">" << "seq2";
	for (unsigned long int i=0; i<seq2.size(); i++){
		if((i%100)==0)
			fastafile << "\n";
		fastafile << seq2[i];
	}
	fastafile <<"\n";

	fastafile.close();
}

//reads the genomic sequence from a file
vector<char> readFASTASeqFromFile (string fileName)
{
	vector<char> seq;
	char c='0';
	ifstream inputF(fileName.c_str());
	if (inputF.is_open())
	  {
		string firstLine;
		getline(inputF, firstLine);
	    while (inputF.get(c))
	    {
	      if(c!='\n')
	    	  seq.push_back(c);
	    }
	    inputF.close();
	  }

	  else cout << "Unable to open file";

	return seq;

}


vector<char> indel (vector<char> ancestorSeq, double indelRate)
{
	//ds: descending sequence
	vector<char> ds = ancestorSeq;

	//ratio of insertion vs deletion
	const double ratio = 1.0;


	random_device rd1;
	mt19937 mt1(rd1());

	//generate random double no to determine insertion/deletion
	random_device rd2;
	mt19937 mt2(rd2());
	uniform_real_distribution<double> dist2(0.0, ratio + 1.0);

	//random number generation to find the substituting base if insertion
	random_device rd3;
	mt19937 mt3(rd3());
	uniform_int_distribution<unsigned long int> dist3(0, 3);

	//generate random int no to determine how long insertion/deletion is (here between 1 and 100)
	random_device rd4;
	mt19937 mt4(rd4());
	uniform_int_distribution<unsigned long int> dist4(1, 100);


	unsigned long int index = 0;
	double change = 0;
	char x = 0;
	int length = 0;

    vector<char> bases = {'A', 'T', 'C', 'G'};

	for(unsigned long int i=0 ; i < indelRate* ancestorSeq.size() ; i++)
	{
		//random number generation to find the site of indel

		uniform_int_distribution<unsigned long int> dist1(0, ds.size()-1);


		index = dist1(mt1);		//indel site
		change = dist2(mt2);
		length = dist4(mt4);

		//in/del based on the probabilities indicated by ratio
		if (change < ratio) //insertion
		{
			//int realIndex = convertIndex(index, alignment[0]);
			vector<char> insertion;
			insertion.reserve(length);
			for(int j=0 ; j<length ; j++)
			{
				insertion.push_back(bases[dist3(mt3)]);
				//alignment[0].insert(alignment[0].begin() + realIndex, x);
				//alignment[1].insert(alignment[1].begin() + realIndex, '-');

			}

			ds.insert(ds.begin() + index, insertion.begin(), insertion.end());

		}
		else {
			if ((change >= ratio)){  //deletion
				ds.erase(ds.begin()+ index , ds.begin() + min(index + length, ds.size()) );

				//int realIndex = convertIndex(index, alignment[0]);
				//for(int j=0; j<length; j++)
				//	alignment[0][realIndex+j]='-';

			}
		}

		//cout<<i<<endl;
	}

	//cout << "\n";

	return ds;
}

int convertIndex (int index,  vector<char> seq)
{
	int counter=0;
	for(unsigned int i=0; i<seq.size(); i++){
		if (seq[i]=='-'){
			continue;
		}

		if (counter==index)
			return i;
		counter++;
	}
	return 0;
}

void printAlignment (vector<vector<char>> alignment, string fileName)
{
	ofstream alignmentfile;
	alignmentfile.open (fileName.c_str());

	for (unsigned long int i=0; i<(alignment[1].size()/100); i++){
		for (int j=0; j<100; j++)
			alignmentfile << alignment[0][i*100 + j];
		alignmentfile << "\n";
		for (int j=0; j<100; j++)
			alignmentfile << alignment[1][i*100 + j];

		alignmentfile << "\n\n";
	}

	for (unsigned int j=0; j<(alignment[1].size()%100); j++)
		alignmentfile << alignment[0][(alignment[1].size()/100)*100 + j];
	alignmentfile << "\n";
	for (unsigned int j=0; j<(alignment[1].size()%100); j++)
			alignmentfile << alignment[1][(alignment[1].size()/100)*100 + j];

	alignmentfile << "\n";


	alignmentfile.close();

}


vector<char> cutSequence (vector<char> seq, double startToCut, double cutRatio)
{
	int start = seq.size()*startToCut;
	int length = seq.size()*cutRatio;
	vector<char> seq2(seq.begin() + start, seq.begin() + start + length -1);
	return seq2;
}


vector<char> repeat (vector<char> ancestorSeq, double repeatRate)
{
	//ds: descending sequence
	vector<char> ds = ancestorSeq;

	//find the origin of repeat
	random_device rd1;
	mt19937 mt1(rd1());
	uniform_int_distribution<unsigned long int> dist1(0, ds.size()-1001);

	//find the site of insertion of the repeat
	random_device rd2;
	mt19937 mt2(rd2());
	uniform_int_distribution<unsigned long int> dist2(0, ds.size()-1);

	//generate random int no to determine how long the repeat is
	random_device rd3;
	mt19937 mt3(rd3());
	uniform_int_distribution<unsigned int> dist3(100, 1000);




	unsigned long int index1 = 0;
	unsigned long int index2 = 0;

	unsigned int length = 0;

    for(unsigned long int i=0 ; i < repeatRate* ancestorSeq.size() ; i++)
	{

		index1 = dist1(mt1);		//repeat site
		index2 = dist2(mt2);		//insertion site
		length = dist3(mt3);		//repeat length

		ds.insert(ds.begin() + index2, ds.begin() + index1, ds.begin() + index1 + length );

	}


	return ds;
}

