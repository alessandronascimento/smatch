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

using namespace std;

class Parser {
public:
	string param;
	string reference_file;
	vector<int> selected_residues;
	vector<string> residue_types;
	vector<string> chain;
	double search_radius;
	string multifile;
	string output_prefix;
	bool write_pdb;
    int matching_residues;
	Parser(char* arg);
	virtual ~Parser();
	void parse_param(string param, ifstream &input);
};

#endif /* PARSER_H_ */
