/*
 * Parser.h
 *
 *  Created on: Oct 25, 2013
 *      Author: asn
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>

using namespace std;

class Parser {
public:
	string param;
	string reference_file;
	vector<int> selected_residues;
	vector<string> residue_types;
	vector<string> chain;
    vector<string> lookup;
	double search_radius;
	string multifile;
	string output_prefix;
    string directory;
	bool write_pdb;
	bool verbose;
    int matching_residues;
    int parallel_jobs;
	Parser(char* arg);
	virtual ~Parser();
	void parse_param(string param, ifstream &input);
};

#endif /* PARSER_H_ */
