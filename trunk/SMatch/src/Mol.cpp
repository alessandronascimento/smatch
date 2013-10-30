/*
 * Mol.cpp
 *
 *  Created on: Sep 30, 2013
 *      Author: asn
 */

#include "Mol.h"

Mol::Mol() {
}

bool Mol::read_pdb(string pdbin){
#ifdef DEBUG
	printf("Opening file %s...\n", pdbin.c_str());
#endif
	if ((pdbin.substr(pdbin.size()-3, 3) == ".gz") or (pdbin.substr(pdbin.size()-2, 2) == ".z")){
		this->read_gzpdb(pdbin);
	}
	else {

		this->filename = pdbin;

		FILE* pdbfile;
		char atom[10];
		char atomname[10];
		char resname[10];
		char chain[2];
		char at[2];
		int itmp;
		int resnumber;
		vector<double> ixyz;
		float  occ, bf, x, y, z;
		pdbfile = fopen(pdbin.c_str(), "r");

		if (pdbfile == NULL){
			printf("Could not open PDB file %s. Please check.\n", pdbin.c_str());
			exit(1);
		}

		while (str[0] != 'A' or str[1] != 'T' or str[2] != 'O' or str[3] != 'M'){ 	//ignoring header
			fgets(str, 80, pdbfile);
		}
		int old_res;

		sscanf(str, "%6s%5d%5s%3s%1s%4d %8f%8f%8f%6f%6f%s", atom, &itmp, atomname, resname, chain, &resnumber, &x, &y, &z, &occ, &bf, at);

		old_res=resnumber;

		Residue* Res = new Residue;
		Res->Natoms=0;

		while (!feof(pdbfile)){
			if(resnumber == old_res and string(atom) == "ATOM"){
				ixyz.push_back(x);
				ixyz.push_back(y);
				ixyz.push_back(z);
				Res->xyz.push_back(ixyz);
				ixyz.clear();
				Res->atomnames.push_back(string(atomname));
				this->check_side_chain(Res, string(atomname));
				Res->Natoms++;
				Res->resname=string(resname);						//redundant
				this->check_restype(Res, string(resname));			//redundant
				Res->resnumber = resnumber;
				Res->chain = string(chain);
//				printf("%4s %4d %4s %6d\n", atomname, itmp, resname, resnumber);
			}
			if (resnumber > old_res and string(atom) == "ATOM"){
				mymol.push_back(*Res);
				delete Res;
				Res = new Residue;

				old_res=resnumber;

				ixyz.push_back(x);
				ixyz.push_back(y);
				ixyz.push_back(z);
				Res->xyz.push_back(ixyz);
				ixyz.clear();
				Res->atomnames.push_back(string(atomname));
				this->check_side_chain(Res, string(atomname));
				Res->Natoms++;
				Res->resname=string(resname);
				this->check_restype(Res, string(resname));
				Res->resnumber = resnumber;
				Res->chain = string(chain);
//				printf("%4s %4d %4s %6d\n", atomname, itmp, resname, resnumber);

			}
				fscanf(pdbfile, "%6s%5d%4s%4s%1s%4d %8f%8f%8f%6f%6f%s", atom, &itmp, atomname, resname, chain, &resnumber, &x, &y, &z, &occ, &bf, at);
			}
		mymol.push_back(*Res);
		fclose(pdbfile);

#ifdef DEBUG
		printf("The molecule has %d residues\n", int(mymol.size()));
		for (unsigned j=0; j< mymol.size(); j++){
			for (unsigned i=0; i< mymol[j].atomnames.size(); i++){
				printf("%4s %4d %4s %4d %8.3f %8.3f %8.3f\n", mymol[j].atomnames[i].c_str(), i+1, mymol[j].resname.c_str(), j+1, mymol[j].xyz[i][0], mymol[j].xyz[i][1],
					mymol[j].xyz[i][2]);
			}
		}
#endif

	}
    this->copy_coordinates();
	return true;
}

bool Mol::read_gzpdb(string pdbin){
	FILE* pdbfile;
	pdbfile = popen(("zcat " + pdbin).c_str(), "r");

	this->filename = pdbin;
	char atom[10];
	char atomname[10];
	char resname[10];
	char chain[2];
	char at[2];
	int itmp;
	int resnumber;
	vector<double> ixyz;
	float  occ, bf, x, y, z;
	string atomtmp;

	if (pdbfile == NULL){
		printf("Could not open PDB file %s. Please check.\n", pdbin.c_str());
		exit(1);
	}

	while (str[0] != 'A' or str[1] != 'T' or str[2] != 'O' or str[3] != 'M'){ 	//ignoring header
		fgets(str, 80, pdbfile);
	}

	int old_res;

	sscanf(str, "%6s%5d%5s%3s%1s%4d %8f%8f%8f%6f%6f%s", atom, &itmp, atomname, resname, chain, &resnumber, &x, &y, &z, &occ, &bf, at);

	old_res=resnumber;

	Residue* Res = new Residue;
	Res->Natoms=0;

	while (!feof(pdbfile)){
		if(resnumber == old_res and string(atom) == "ATOM"){
			ixyz.push_back(x);
			ixyz.push_back(y);
			ixyz.push_back(z);
			Res->xyz.push_back(ixyz);
			ixyz.clear();
			Res->atomnames.push_back(string(atomname));
			this->check_side_chain(Res, string(atomname));
			Res->Natoms++;
			Res->resname=string(resname);						//redundant
			this->check_restype(Res, string(resname));			//redundant
			Res->resnumber = resnumber;
			Res->chain = string(chain);
			//				printf("%4s %4d %4s %6d\n", atomname, itmp, resname, resnumber);
		}
		if (resnumber > old_res and string(atom) == "ATOM"){
			mymol.push_back(*Res);
			delete Res;
			Res = new Residue;

			old_res=resnumber;

			ixyz.push_back(x);
			ixyz.push_back(y);
			ixyz.push_back(z);
			Res->xyz.push_back(ixyz);
			ixyz.clear();
			Res->atomnames.push_back(string(atomname));
			this->check_side_chain(Res, string(atomname));
			Res->Natoms++;
			Res->resname=string(resname);
			this->check_restype(Res, string(resname));
			Res->resnumber = resnumber;
			Res->chain = string(chain);
			//				printf("%4s %4d %4s %6d\n", atomname, itmp, resname, resnumber);

		}
		fscanf(pdbfile, "%6s%5d%4s%4s%1s%4d %8f%8f%8f%6f%6f%s", atom, &itmp, atomname, resname, chain, &resnumber, &x, &y, &z, &occ, &bf, at);
	}
	mymol.push_back(*Res);

#ifdef DEBUG
	printf("The molecule has %d residues\n", int(mymol.size()));
	for (unsigned j=0; j< mymol.size(); j++){
		for (unsigned i=0; i< mymol[j].atomnames.size(); i++){
			printf("%4s %4d %4s %4d %8.3f %8.3f %8.3f\n", mymol[j].atomnames[i].c_str(), i+1, mymol[j].resname.c_str(), j+1, mymol[j].xyz[i][0], mymol[j].xyz[i][1],
					mymol[j].xyz[i][2]);
		}
	}
#endif

	pclose(pdbfile);
    this->copy_coordinates();
	return true;
}

Mol::~Mol() {
	// TODO Auto-generated destructor stub
}

char Mol::residues_3_to_1_letter(string residue){
	char res = 'X';
	if (residue == "ALA"){
		res = 'A';
	}
	else if (residue == "ARG"){
		res = 'R';
	}
	else if (residue == "ASN"){
		res = 'N';
	}
	else if (residue == "ASP"){
		res = 'D';
	}
	else if (residue == "CYS"){
		res = 'C';
	}
	else if (residue == "GLU"){
		res = 'E';
	}
	else if (residue == "GLN"){
		res = 'Q';
	}
	else if (residue == "GLY"){
		res = 'G';
	}
	else if (residue == "HIS"){
		res = 'H';
	}
	else if (residue == "ILE"){
		res = 'I';
	}
	else if (residue == "LEU"){
		res = 'L';
	}
	else if (residue == "LYS"){
		res = 'K';
	}
	else if (residue == "MET"){
		res = 'M';
	}
	else if (residue == "PHE"){
		res = 'F';
	}
	else if (residue == "PRO"){
		res = 'P';
	}
	else if (residue == "SER"){
		res = 'S';
	}
	else if (residue == "THR"){
		res = 'T';
	}
	else if (residue == "TRP"){
		res = 'W';
	}
	else if (residue == "TYR"){
		res = 'Y';
	}
	else if (residue == "VAL"){
		res = 'V';
	}
	else {
		res = 'X';
	}
	return res;
}

void Mol::check_side_chain(Residue* Res, string atomname){
	if (atomname == "CB"){
		Res->has_CB = true;
	}
	else if (atomname == "CG"){
		Res->has_CG = true;
	}
	else if (atomname == "CD"){
		Res->has_CD = true;
	}
}

void Mol::check_restype(Residue* Res, string resname){
	if (resname == "GLU" or resname == "ASP"){
		Res->is_acid = true;
	}
	else {
		Res->is_acid = false;
	}
	if (resname == "LYS" or resname == "ARG"){
		Res->is_basic = true;
	}
	else {
		Res->is_basic = false;
	}
	if (resname == "ALA" or resname == "GLY" or resname == "VAL" or resname == "LEU" or resname == "ILE" or resname == "PHE"){
		Res->is_apolar = true;
	}
	else {
		Res->is_apolar = false;
	}
}

void Mol::copy_coordinates(void){
    for (unsigned i=0; i<this->mymol.size(); i++){
        for (unsigned j=0; j< this->mymol[i].xyz.size(); j++){
            this->xyz.push_back(this->mymol[i].xyz[j]);
        }
    }
}
