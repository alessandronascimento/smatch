/*
 * Engine.h
 *
 *  Created on: Oct 29, 2013
 *      Author: asn
 */

#ifndef ENGINE_H_
#define ENGINE_H_

#include <stdio.h>
#include <stdlib.h>
#include "Mol.h"
#include "Optimization.h"
#include "SMatch.h"
#include "Printer.h"
#include "Parser.h"
#include "SMatch.h"
#ifdef HAS_MPI
#include "boost/mpi.hpp"
namespace mpi = boost::mpi;
#endif

using namespace std;
class Engine {
public:
    Printer* Writer;
    Engine(Printer* _Writer);
	virtual ~Engine();
	void mol_extraction(Mol* M1, Mol* MExtract1, Parser* Input);
	void mol_extraction(Mol* M1, Mol* MExtract1, string resname);
	vector<string> make_unique(vector<string> resnames);
    int serial_search(Mol* ME1, Parser* Input, vector<string> unique, vector<string> pdb_list);
    int run_serial(Mol* ME1, vector<string> unique, Parser* Input);
    int serial_search_omp(Mol* ME1, Parser* Input, vector<string> unique, vector<string> pdb_list);
	void print_mol_info(Mol* M);
	bool residues_match(string r1, string r2);
	void mol_extraction(Mol* M1, Mol* MExtract1, vector<string> resnames);
    char info[98];
#ifdef HAS_MPI
    int run_over_mpi(int argc, char* argv[], Mol* ME1, vector<string> unique, Parser* Input);
#endif
};

#endif /* ENGINE_H_ */
