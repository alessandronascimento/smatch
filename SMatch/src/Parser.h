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

using namespace std;

class Parser {
	string referece_file;
	vector<int> selected_residues;
	vector<string> residue_types;

public:
	Parser();
	virtual ~Parser();
};

#endif /* PARSER_H_ */
