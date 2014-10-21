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
	verbose = false;
    matching_residues = -1;
    this->directory = "./";
	ifstream input(arg);
	char line[256];
    if (! input.is_open()){
        printf("Could not open file %s. Please check.\n", arg);
        exit(1);
    }
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
    string lu;
	if (param == "reference_file"){
		input >> this->reference_file;
	}
	else if (param == "multifile"){
		input >> this->multifile;
	}
	else if (param == "search"){
        input >> res >> ch >> resname >> lu;
		this->selected_residues.push_back(res);
		this->residue_types.push_back(resname);
		this->chain.push_back(ch);
        this->lookup.push_back(lu);
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
    else if (param == "directory"){
        input >> this->directory;
    }
    else if (param == "verbose"){
    	input >> lu;
    	if ((lu == "yes") or (lu == "YES") or (lu == "Yes")){
    		this->verbose = true;
    	}
    }
    else if ( param == "parallel_jobs"){
        input >> this->parallel_jobs;
    }
}

Parser::~Parser() {
}
