/*
 * Engine.cpp
 *
 *  Created on: Oct 29, 2013
 *      Author: asn
 */

#include "Engine.h"

Engine::Engine() {
}

void Engine::print_mol_info(Mol* M){
	printf("Mol %12.12s has %d residues matching the criteria:\n", M->filename.c_str(), int(M->mymol.size()));
	for (unsigned i=0; i< M->mymol.size(); i++){
		printf("\t %3.3s %1.1s %4d\n", M->mymol[i].resname.c_str(), M->mymol[i].chain.c_str(), M->mymol[i].resnumber);
	}
}

Engine::~Engine() {
	printf("Quiting SMatch...\n");
}

void Engine::mol_extraction(Mol* M1, Mol* MExtract1, Parser* Input){
	for (unsigned i=0; i<Input->selected_residues.size(); i++){
		for (unsigned j=0; j< M1->mymol.size(); j++){
			if ((M1->mymol[j].resnumber == Input->selected_residues[i]) and (M1->mymol[j].chain == Input->chain[i])){
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
            if ((resnames[i] == unique[j]) or (resnames[i]=="*")){
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
	for (unsigned i=0; i< pdb_list.size(); i++) {
        Mol* M2 = new Mol(Input);
        Mol* MExtract2 = new Mol(Input);
		if (M2->read_pdb(pdb_list[i])){
			for (unsigned j=0; j< unique.size(); j++){
				mol_extraction(M2, MExtract2, unique[j]);
			}

#ifdef DEBUG
			printf("Found %d matching residues in search molecule.\n", int(MExtract2->mymol.size()));
#endif

			vector<vector<vector<double> > >xyz;
            Optimization* Opt = new Optimization(Writer, ME1, Input);
			opt_result_t* opt_result = new opt_result_t;
			Opt->optimize_rmsd(MExtract2, opt_result);

            if (opt_result->succeded){
                printf("\tMatched residues:\n");
                for (unsigned a=0; a< opt_result->imatched1.size(); a++){
                    printf("\t\t %15.15s %4.4s %4d       %15.15s %4.4s %4d %6.4f\n", ME1->filename.c_str(), opt_result->smatched1[a].c_str(), opt_result->imatched1[a],
                           MExtract2->filename.c_str(), opt_result->smatched2[a].c_str(), opt_result->imatched2[a], opt_result->rmsds[a]);
                }
            }

            else {
                printf("No matched results for %-15.15s\n", MExtract2->filename.c_str());
            }

			if (opt_result->succeded and Input->write_pdb){
				Coord* CoordManip = new Coord;
				xyz = CoordManip->rototranslate(M2, opt_result->rotrans[0], opt_result->rotrans[1], opt_result->rotrans[2], opt_result->rotrans[3],
						opt_result->rotrans[4], opt_result->rotrans[5]);
				delete CoordManip;
				Writer->write_pdb(M2, xyz, 0.0, opt_result->rmsd, (pdb_list[i].substr(0,pdb_list[i].find(".pdb")) + "_smatch"));
			}
			delete Opt;
			delete opt_result;
		}
		delete MExtract2;
		delete M2;
	}
	return EXIT_SUCCESS;
}

#ifdef HAS_MPI
int Engine::run_over_mpi(int argc, char* argv[], Mol* ME1, vector<string> unique, Parser* Input, Printer* Writer){
	mpi::environment env(argc, argv);
	mpi::communicator world;
	vector<string> pdb_list;
	vector<string> tmp;
	vector<vector<string> > chuncks(world.size());

	if (world.rank() == 0 ){
        Writer->print_welcome();
        this->print_mol_info(ME1);

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

//		printf("There are %d files, divided into %d chuncks of size %d.\n", int(pdb_list.size()), world.size(), chunck_size);

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

#ifdef DEBUG
		for (unsigned i=0; i< chuncks.size(); i++){
			printf("Chuncks[%d] has %d elements.\n", int(i), int(chuncks[i].size()));
		}
#endif

		scatter(world, chuncks, tmp,0);
		serial_search(ME1, Writer, Input, unique, tmp);
	}
	else{
		scatter(world, tmp, 0);
		serial_search(ME1, Writer, Input, unique, tmp);
	}
	return 0;
}
#endif

int Engine::run_serial(Mol* ME1, vector<string> unique, Parser* Input, Printer* Writer){
	vector<string> pdb_list;
	ifstream multifile;
	multifile.open(Input->multifile.c_str());
	string pdbmol;
	multifile >> pdbmol;

	while (!multifile.eof() and pdbmol != "EOF"){
        pdb_list.push_back(pdbmol);
		multifile >> pdbmol;
	}
	multifile.close();

	this->serial_search_omp(ME1, Writer, Input, unique, pdb_list);

	return 0;

}

int Engine::serial_search_omp(Mol* ME1, Printer* Writer, Parser* Input, vector<string> unique, vector<string> pdb_list) {

#pragma omp parallel
	{
#pragma omp for schedule(static, 1)
		for (unsigned i=0; i< pdb_list.size(); i++) {
            Mol* M2 = new Mol(Input);
            Mol* MExtract2 = new Mol(Input);
			if (M2->read_pdb(pdb_list[i])){
				for (unsigned j=0; j< unique.size(); j++){
					mol_extraction(M2, MExtract2, unique[j]);
				}

#ifdef DEBUG
				printf("Found %d matching residues in search molecule.\n", int(MExtract2->mymol.size()));
#endif

				vector<vector<vector<double> > >xyz;
                Optimization* Opt = new Optimization(Writer, ME1, Input);
				opt_result_t* opt_result = new opt_result_t;
//				Opt->optimize_rmsd(MExtract2, opt_result);
                Opt->optimize_rmsd(M2, opt_result);

                if (opt_result->succeded){
                    printf("\tMatched residues:\n");
                    for (unsigned a=0; a< opt_result->imatched1.size(); a++){
                        printf("\t\t %15.15s %4.4s %4d       %15.15s %4.4s %4d %6.4f\n", ME1->filename.c_str(), opt_result->smatched1[a].c_str(), opt_result->imatched1[a],
                               MExtract2->filename.c_str(), opt_result->smatched2[a].c_str(), opt_result->imatched2[a], opt_result->rmsds[a]);
                    }
                }
                else {
                    printf("No matched results for %-15.15s\n", MExtract2->filename.c_str());
                }

				if (opt_result->succeded and Input->write_pdb){
					Coord* CoordManip = new Coord;
					xyz = CoordManip->rototranslate(M2, opt_result->rotrans[0], opt_result->rotrans[1], opt_result->rotrans[2], opt_result->rotrans[3],
							opt_result->rotrans[4], opt_result->rotrans[5]);
					delete CoordManip;
					Writer->write_pdb(M2, xyz, 0.0, opt_result->rmsd, (pdb_list[i].substr(0,pdb_list[i].find(".pdb")) + "_smatch"));
				}
				delete Opt;
				delete opt_result;
			}
			delete MExtract2;
			delete M2;
		}
	}
	return EXIT_SUCCESS;
}

vector<string> Engine::check_resnames(vector<string> resnames){
    vector<string> all_residues;
    bool has_all = false;
    bool has_acid = false;
    bool has_basic = false;
    for (unsigned i=0; i<resnames.size(); i++){
        if (resnames[i] == "*"){
            has_all = true;
        }
        else if (resnames[i] == "acid"){
            has_acid = true;
        }
        else if (resnames[i] == "basic"){
            has_basic = true;
        }
    }
    if (has_all){
        all_residues.push_back("ARG"); // 1
        all_residues.push_back("LYS"); // 2
        all_residues.push_back("HIS"); // 3
        all_residues.push_back("ASP"); // 4
        all_residues.push_back("GLU"); // 5
        all_residues.push_back("SER"); // 6
        all_residues.push_back("THR"); // 7
        all_residues.push_back("ASN"); // 8
        all_residues.push_back("GLN"); // 9
        all_residues.push_back("CYS"); // 10
        all_residues.push_back("GLY");
        all_residues.push_back("PRO");
        all_residues.push_back("ALA");
        all_residues.push_back("VAL");
        all_residues.push_back("ILE"); // 15
        all_residues.push_back("LEU");
        all_residues.push_back("MET");
        all_residues.push_back("PHE");
        all_residues.push_back("TYR");
        all_residues.push_back("TRP"); // 20
    }
    if (has_acid){
        all_residues.push_back("ASP"); // 4
        all_residues.push_back("GLU");
    }
    if (has_basic){
        all_residues.push_back("ARG"); // 1
        all_residues.push_back("LYS");
    }

    return (all_residues);
}
