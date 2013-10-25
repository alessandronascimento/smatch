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

Mol* Optimization::M1;

void mol_extraction(Mol* M1, Mol* MExtract1, char* argv[], int number_of_selected_res){
	int resn;
	for (int i=0; i< number_of_selected_res; i++){
		resn = atoi(argv[2+i]);
		resn = resn-M1->resnumber_phase;
		MExtract1->resnames.push_back(M1->resnames[resn]);
		MExtract1->residue_pointer.push_back(int(MExtract1->xyz.size()));
		printf("Residue %d selected (%s): %dth residue\n", atoi(argv[2+i]), M1->resnames[resn].c_str(), resn);
			for (int j=M1->residue_pointer[resn]; j<M1->residue_pointer[resn+1]; j++){
				MExtract1->xyz.push_back(M1->xyz[j]);
				MExtract1->atomnames.push_back(M1->atomnames[j]);
			}
		}
	MExtract1->N = int(MExtract1->xyz.size());
	MExtract1->Nres = int(MExtract1->resnames.size());
	MExtract1->filename = M1->filename;
//	MExtract1->charge_distribution();

#ifdef DEBUG
	for (unsigned i=0; i<MExtract1->resnames.size()-1; i++){
		printf("%s\n", MExtract1->resnames[i].c_str());
		for (int j=MExtract1->residue_pointer[i]; j<MExtract1->residue_pointer[i+1]; j++){
			printf("Atom %5s   %8.3f %8.3f %8.3f    %8.3f\n", MExtract1->atomnames[j].c_str(), MExtract1->xyz[j][0], MExtract1->xyz[j][1], MExtract1->xyz[j][2], MExtract1->charges[j]);
		}
	}
	printf("%s\n", MExtract1->resnames[MExtract1->resnames.size()-1].c_str());
	int j=MExtract1->residue_pointer[MExtract1->resnames.size()-1];
	while (j<MExtract1->N){
		printf("Atom %5s   %8.3f %8.3f %8.3f    %8.3f\n", MExtract1->atomnames[j].c_str(), MExtract1->xyz[j][0], MExtract1->xyz[j][1], MExtract1->xyz[j][2], MExtract1->charges[j]);;
		j++;
	}
#endif
}

void mol_extraction(Mol* M1, Mol* MExtract1, string resname){
	int start_atom, end_atom;
	int rpointer=0;
	for (unsigned i=0; i< M1->resnames.size(); i++){
		if (M1->resnames[i] == resname){

			MExtract1->resnames.push_back(M1->resnames[i]);
			start_atom = M1->residue_pointer[i];
			if (i < M1->resnames.size()-1){
				end_atom = M1->residue_pointer[i+1];
			}
			else{
				end_atom = M1->N;
			}
			
			MExtract1->residue_pointer.push_back(rpointer);
			for (int j=start_atom; j<end_atom; j++){
				MExtract1->xyz.push_back(M1->xyz[j]);
				MExtract1->atomnames.push_back(M1->atomnames[j]);
				rpointer++;
			}
//#ifdef DEBUG
			printf("%3d -  residue %s [%d]\n", MExtract1->resnames.size()-1, M1->resnames[i].c_str(), i+M1->resnumber_phase);
//#endif
		}
	}

	MExtract1->N = int(MExtract1->xyz.size());
	MExtract1->Nres = int(MExtract1->resnames.size());
	MExtract1->filename = M1->filename;
//	MExtract1->charge_distribution();
}

int main(int argc, char* argv[]) {

	if (argc < 2){
		printf ("Usage: %s structure.pdb res1 res2 res3.\n", argv[0]);
		exit(1);
	}

	Mol* M1 = new Mol;
	M1->read_pdb(string (argv[1]));

	int number_of_selected_res= argc-2;
	printf("%d residues were selected for structure extraction...\n", number_of_selected_res);

	Mol* MExtract1 = new Mol;
	mol_extraction(M1, MExtract1, argv, number_of_selected_res);

	Coord* CoordManip = new Coord;

	vector<string> res_types;
	bool not_found = false;
	res_types.push_back(MExtract1->resnames[0]);
	for (unsigned i =0; i< MExtract1->resnames.size(); i++){
		not_found = false;
		for (unsigned j=0; j< res_types.size(); j++){
			if (MExtract1->resnames[i] == res_types[j]){
				not_found = true;
			}
		}
		if (! not_found){
			res_types.push_back(MExtract1->resnames[i]);
		}
	}

#ifdef DEBUG
	printf("List of residues to look for:\n");
	for (unsigned i=0; i< res_types.size(); i++){
		printf("\t%s\n", res_types[i].c_str());
	}
#endif

	Printer* Writer = new Printer;
	Writer->write_pdb(MExtract1, MExtract1->xyz, 0.0, 0.0, "ME1");

	ifstream multifile;
	multifile.open("INSMATCH");
	string pdbmol;
	multifile >> pdbmol;

	while (!multifile.eof()){
		Mol* M2 = new Mol;
		Mol* MExtract2 = new Mol;
		M2->read_pdb(pdbmol);
		for (unsigned i=0; i< res_types.size(); i++){
			mol_extraction(M2, MExtract2, res_types[i]);
		}

		printf("Found %d matching residues in search molecule.\n", int(MExtract2->resnames.size()));

		vector<vector<double> > xyz;
		Optimization* Opt = new Optimization(Writer, MExtract1);
		opt_result_t* opt_result = new opt_result_t;
		Opt->optimize_rmsd(MExtract2, opt_result);

		if (opt_result->succeded){
			M2->xyz = CoordManip->rototranslate(M2->xyz, M2, opt_result->rotrans[0], opt_result->rotrans[1], opt_result->rotrans[2], opt_result->rotrans[3],
						opt_result->rotrans[4], opt_result->rotrans[5]);
			Writer->write_pdb(M2, M2->xyz, 0.0, opt_result->rmsd, pdbmol);
			Writer->write_pdb(MExtract2, opt_result->xyz, 0.0, 0.0, "ME2");
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
