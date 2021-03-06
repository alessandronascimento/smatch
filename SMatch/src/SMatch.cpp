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

    char info[98];
	Parser* Input = new Parser(argv[1]);
    Printer* Writer = new Printer(Input);

    Engine* RunEngine = new Engine(Writer);
    Mol* M1 = new Mol(Input);
    sprintf(info, "Reading reference file %s...", Input->reference_file.c_str());
    Writer->print_info(info);
    M1->read_pdb(Input->reference_file);

    Mol* MExtract1 = new Mol(Input);
	RunEngine->mol_extraction(M1, MExtract1, Input);
    Writer->write_pdb(MExtract1, 0.0, 0.0, "ME1");
	vector<string> unique = RunEngine->make_unique(Input->lookup);

    sprintf(info, "Searching for matches...");
    Writer->print_info(info);

#ifdef HAS_MPI
    RunEngine->run_over_mpi(argc, argv, MExtract1, unique, Input);
#else
    RunEngine->run_serial(MExtract1, unique, Input);
#endif

	unique.clear();

	delete RunEngine;
	delete MExtract1;
	delete M1;
	delete Input;
	delete Writer;

	return EXIT_SUCCESS;
}
