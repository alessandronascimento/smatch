//============================================================================
// Name        : SMatch.cpp
// Author      : Alessandro Nascimento
// Version     :
// Copyright   : Rights Reserved
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "Mol.h"
#include "Optimization.h"
#include "SMatch.h"
#include "Printer.h"
#include "Parser.h"

Mol* Optimization::M1;

void mol_extraction(Mol* M1, Mol* MExtract1, Parser* Input){
	for (unsigned i=0; i<Input->selected_residues.size(); i++){
		for (unsigned j=0; j< M1->mymol.size(); j++){
			if (M1->mymol[j].resnumber == Input->selected_residues[i]){
				MExtract1->mymol.push_back(M1->mymol[j]);
			}
		}
	}
}

void mol_extraction(Mol* M1, Mol* MExtract1, string resname){
	for (unsigned i=0; i<M1->mymol.size(); i++){
		if (M1->mymol[i].resname == resname){
			MExtract1->mymol.push_back(M1->mymol[i]);
		}
	}
}

vector<string> make_unique(vector<string> resnames){
	vector<string> unique;
	unique.push_back(resnames[0]);
	bool found;
	for (unsigned i=1; i<resnames.size(); i++){
		found=false;
		for (unsigned j=0; j< unique.size(); j++){
			if (resnames[i] == unique[j]){
				found = true;
			}
		}
		if (!found){
			unique.push_back(resnames[i]);
		}
	}
	return(unique);
}

int main(int argc, char* argv[]) {

	if (argc != 2){
		printf ("Usage: %s input_file.\n", argv[0]);
		exit(1);
	}

	Parser* Input = new Parser(argv[1]);
    Coord* CoordManip = new Coord;

	Mol* M1 = new Mol;
	M1->read_pdb(Input->reference_file);

	int number_of_selected_res= int(Input->selected_residues.size());
	printf("%d residues were selected for structure extraction...\n", number_of_selected_res);

	Mol* MExtract1 = new Mol;
	mol_extraction(M1, MExtract1, Input);

	vector<string> unique = make_unique(Input->residue_types);

	Printer* Writer = new Printer;
	Writer->write_pdb(MExtract1, 0.0, 0.0, "ME1");

	ifstream multifile;
	multifile.open(Input->multifile.c_str());
	string pdbmol;
	multifile >> pdbmol;

	while (!multifile.eof()){
		Mol* M2 = new Mol;
		Mol* MExtract2 = new Mol;
		M2->read_pdb(pdbmol);
		for (unsigned i=0; i< unique.size(); i++){
			mol_extraction(M2, MExtract2, unique[i]);
		}

		printf("Found %d matching residues in search molecule.\n", int(MExtract2->mymol.size()));

        vector<vector<vector<double> > >xyz;
        Optimization* Opt = new Optimization(Writer, MExtract1);
		opt_result_t* opt_result = new opt_result_t;
		Opt->optimize_rmsd(MExtract2, opt_result);

		if (opt_result->succeded){
            xyz = CoordManip->rototranslate(M2, opt_result->rotrans[0], opt_result->rotrans[1], opt_result->rotrans[2], opt_result->rotrans[3],
						opt_result->rotrans[4], opt_result->rotrans[5]);
//			Writer->write_pdb(M2, M2->xyz, 0.0, opt_result->rmsd, pdbmol);
//			Writer->write_pdb(MExtract2, opt_result->xyz, 0.0, 0.0, "ME2");
		}
		delete MExtract2;
		delete M2;
		multifile >> pdbmol;
	}


	delete Writer;
	delete MExtract1;
	delete M1;
    return EXIT_SUCCESS;
}





//TODO: Alter SMatch.cpp for the new definition of a molecule in Mol class.
