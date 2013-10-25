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

using namespace std;

struct Residue{
	vector<vector<double> > xyz;
	vector<string> atomnames;
	string resname;
	char chain;
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

	void check_side_chain(Residue* Res, string atomname);
	void check_restype(Residue* Res, string resname);
	char str[80];
	string filename;
	Mol();
	bool read_pdb(string pdbin);
	bool read_gzpdb(string pdbin);
	char residues_3_to_1_letter(string residue);
	virtual ~Mol();
};

#endif /* MOL_H_ */
