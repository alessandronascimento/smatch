//============================================================================
// Name        : SMatch.cpp
// Author      : Alessandro Nascimento
// Version     :
// Copyright   : Rights Reserved
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "Mol.h"
#include "Optimization.h"
#include "Printer.h"
#include "Parser.h"
#include "Engine.h"

Mol* Optimization::M1;

int main(int argc, char* argv[]) {

	if (argc != 2){
		printf ("Usage: %s input_file.\n", argv[0]);
		exit(1);
	}

	Parser* Input = new Parser(argv[1]);

	Engine* RunEngine = new Engine;

	Mol* M1 = new Mol;
	M1->read_pdb(Input->reference_file);

	int number_of_selected_res= int(Input->selected_residues.size());
	printf("%d residues were selected for structure extraction...\n", number_of_selected_res);

	Mol* MExtract1 = new Mol;
	RunEngine->mol_extraction(M1, MExtract1, Input);


	vector<string> unique = RunEngine->make_unique(Input->residue_types);

	Printer* Writer = new Printer(Input);
	Writer->write_pdb(MExtract1, 0.0, 0.0, "ME1");

	RunEngine->run_over_mpi(argc, argv, MExtract1, unique, Input, Writer);

	delete RunEngine;
	delete MExtract1;
	delete M1;
	delete Input;
	delete Writer;

	return EXIT_SUCCESS;
}
