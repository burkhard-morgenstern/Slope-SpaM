/*
 * Header.h
 *
 *  Created on: 2 May 2016
 *      Author: salma
 */
//header file for generator.cpp
//Sequence evolution simulation
#ifndef HEADER_H_
#define HEADER_H_


#include <iostream>
//#include <iomanip>
//#include <string>
//#include <map>
#include <random>
#include <bits/random.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <time.h>
 #include <getopt.h>

using namespace std;

vector<char> DNA_SeqGenerator(unsigned long int);
vector<char> mutate(vector<char>, double);
vector<char> mutate2(vector<char>, double);
double Jukes_Cantor(vector<char>&, vector<char>&);
void Print (vector<char>&);
void saveFASTA (vector<char>, vector<char>, string);
vector<char> readFASTASeqFromFile (string);
vector<char> indel (vector<char> , double);
int convertIndex (int ,  vector<char>);
void printAlignment (vector<vector<char>> , string);
vector<char> cutSequence (vector<char> , double , double);
vector<char> repeat (vector<char>, double);

#endif /* HEADER_H_ */
