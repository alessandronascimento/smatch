/*
 * Engine.cpp
 *
 *  Created on: Oct 29, 2013
 *      Author: asn
 */

#include "Engine.h"

Engine::Engine() {

}

Engine::~Engine() {
}

void Engine::mol_extraction(Mol* M1, Mol* MExtract1, Parser* Input){
	for (unsigned i=0; i<Input->selected_residues.size(); i++){
		for (unsigned j=0; j< M1->mymol.size(); j++){
			if (M1->mymol[j].resnumber == Input->selected_residues[i]){
				MExtract1->mymol.push_back(M1->mymol[j]);
			}
		}
	}
	MExtract1->filename = M1->filename;
}

void Engine::mol_extraction(Mol* M1, Mol* MExtract1, string resname){
	for (unsigned i=0; i<M1->mymol.size(); i++){
		if (M1->mymol[i].resname == resname){
			MExtract1->mymol.push_back(M1->mymol[i]);
		}
	}
	MExtract1->filename = M1->filename;
}

vector<string> Engine::make_unique(vector<string> resnames){
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

int Engine::serial_search(Mol* ME1, Printer* Writer, Parser* Input, vector<string> unique, vector<string> pdb_list) {

	Coord* CoordManip = new Coord;

	for (unsigned i=0; i< pdb_list.size(); i++) {
		Mol* M2 = new Mol;
		Mol* MExtract2 = new Mol;
		M2->read_pdb(pdb_list[i]);
		for (unsigned j=0; j< unique.size(); j++){
			mol_extraction(M2, MExtract2, unique[j]);
		}

#ifdef DEBUG
		printf("Found %d matching residues in search molecule.\n", int(MExtract2->mymol.size()));
#endif

        vector<vector<vector<double> > >xyz;
        Optimization* Opt = new Optimization(Writer, ME1);
		opt_result_t* opt_result = new opt_result_t;
		Opt->optimize_rmsd(MExtract2, opt_result);

		if (opt_result->succeded){
            xyz = CoordManip->rototranslate(M2, opt_result->rotrans[0], opt_result->rotrans[1], opt_result->rotrans[2], opt_result->rotrans[3],
						opt_result->rotrans[4], opt_result->rotrans[5]);
			Writer->write_pdb(M2, xyz, 0.0, opt_result->rmsd, (pdb_list[i].substr(0,pdb_list[i].find(".pdb")) + "_smatch"));
		}
		delete MExtract2;
		delete M2;
	}
    return EXIT_SUCCESS;
}

int Engine::run_over_mpi(int argc, char* argv[], Mol* ME1, vector<string> unique, Parser* Input, Printer* Writer){
	mpi::environment env(argc, argv);
	mpi::communicator world;
	vector<string> pdb_list;
	vector<string> tmp;
	vector<vector<string> > chuncks(world.size());

	if (world.rank() == 0 ){

		ifstream multifile;
		multifile.open(Input->multifile.c_str());
		string pdbmol;
		multifile >> pdbmol;

		while (!multifile.eof() and pdbmol != "EOF"){
			pdb_list.push_back(pdbmol);
			multifile >> pdbmol;
		}
		multifile.close();

		int chunck_size = int(pdb_list.size()/world.size());

		printf("There are %d files, divided into %d chuncks of size %d.\n", int(pdb_list.size()), world.size(), chunck_size);

		for (int i=0; i< world.size(); i++){
			tmp.clear();
			for (int j=0; j< chunck_size; j++){
				pdbmol = pdb_list.front();
				tmp.push_back(pdbmol);
				pdb_list.erase(pdb_list.begin());
			}
			chuncks[i] = tmp;
		}
		tmp.clear();

		for (unsigned i=0; i<pdb_list.size(); i++){
			chuncks[i].push_back(pdb_list[i]);
		}

		scatter(world, chuncks, tmp,0);
		serial_search(ME1, Writer, Input, unique, tmp);
	}
	else{
		scatter(world, tmp, 0);
		serial_search(ME1, Writer, Input, unique, tmp);
	}
	return 0;
}
