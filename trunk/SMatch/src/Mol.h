/*
 * Mol.h
 *
 *  Created on: Sep 30, 2013
 *      Author: asn
 */

#ifndef MOL_H_
#define MOL_H_


#include <fstream>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string.h>
#include <zlib.h>
#include "Parser.h"

using namespace std;

struct Residue{
	vector<vector<double> > xyz;
	vector<string> atomnames;
	string resname;
	string chain;
	int resnumber;
	int Natoms;
	string atom_type;
	bool is_acid;
	bool is_basic;
	bool is_apolar;
	bool has_CB;
	bool has_CG;
	bool has_CD;

};

class Mol {
public:
	vector<Residue> mymol;
    vector<vector<double> >xyz;

	void check_side_chain(Residue* Res, string atomname);
	void check_restype(Residue* Res, string resname);
	char str[100];
	string filename;
    Parser* Input;
    Mol(Parser* _Input);
	bool read_pdb(string pdbin);
	bool read_gzpdb(string pdbin);
	char residues_3_to_1_letter(string residue);
	virtual ~Mol();
    void copy_coordinates(void);
};

#endif /* MOL_H_ */
