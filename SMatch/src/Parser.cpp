/*
 * Parser.cpp
 *
 *  Created on: Oct 25, 2013
 *      Author: asn
 */

#include "Parser.h"


Parser::Parser(char* arg) {
	search_radius=5.0;
	write_pdb = false;
    matching_residues = -1;
	ifstream input(arg);
	char line[256];
	while (!input.eof()){
		input >> param;
		if ((param.substr(0,1) == "#") or (param.substr(0,2) == "//")){
			input.getline(line, 256);
		}
		else{
			this->parse_param(param, input);
		}
	}
	input.close();
}


void Parser::parse_param(string param, ifstream &input){
	int res;
	string resname;
	string ch;
	if (param == "reference_file"){
		input >> this->reference_file;
	}
	else if (param == "multifile"){
		input >> this->multifile;
	}
	else if (param == "search"){
		input >> res >> ch >> resname;
		this->selected_residues.push_back(res);
		this->residue_types.push_back(resname);
		this->chain.push_back(ch);
	}
	else if (param == "search_radius"){
		input >> this->search_radius;
	}
	else if (param == "output_prefix"){
		input >> this->output_prefix;
	}
	else if (param == "write_pdb"){
		input >> ch;
		if (ch == "yes" or ch == "YES" or ch == "Yes"){
			this->write_pdb = true;
		}
	}
    else if (param == "matching_residues"){
        input >> this->matching_residues;
    }
}

Parser::~Parser() {
}
