/*
 * Engine.cpp
 *
 *  Created on: Oct 29, 2013
 *      Author: asn
 */

#include "Engine.h"

Engine::Engine(Printer* _Writer) {
    this->Writer = _Writer;
}


void Engine::print_mol_info(Mol* M){
    sprintf(info, "Mol %12.12s has %d residues matching the criteria:", M->filename.c_str(), int(M->mymol.size()));
    Writer->print_info(info);
	for (unsigned i=0; i< M->mymol.size(); i++){
        sprintf(info, "    %3.3s %1.1s %4d", M->mymol[i].resname.c_str(), M->mymol[i].chain.c_str(), M->mymol[i].resnumber);
        Writer->print_info(info);
	}
}

Engine::~Engine() {
    sprintf(info, "Quiting SMatch...");
    Writer->print_info(info);
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
		if (this->residues_match(M1->mymol[i].resname, resname)){
			MExtract1->mymol.push_back(M1->mymol[i]);
		}
	}
	MExtract1->filename = M1->filename;
}

void Engine::mol_extraction(Mol* M1, Mol* MExtract1, vector<string> resnames){
	bool has_acid=false, has_basic=false, has_polar=false, has_all=false;
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
		else if (resnames[i] == "polar"){
			has_polar = true;
		}
	}

	if (has_acid){
		for (unsigned i=0; i< M1->mymol.size(); i++){
			if (M1->mymol[i].is_acid){
				MExtract1->mymol.push_back(M1->mymol[i]);
			}
		}
	}

	if (has_all){
		for (unsigned i=0; i< M1->mymol.size(); i++){
			MExtract1->mymol.push_back(M1->mymol[i]);
		}
	}

	else if (has_basic){
		for (unsigned i=0; i< M1->mymol.size(); i++){
			if (M1->mymol[i].is_basic){
				MExtract1->mymol.push_back(M1->mymol[i]);
			}
		}
	}

	else if (has_polar){
		for (unsigned i=0; i< M1->mymol.size(); i++){
			if (!M1->mymol[i].is_apolar and !M1->mymol[i].is_acid and !M1->mymol[i].is_basic){
				MExtract1->mymol.push_back(M1->mymol[i]);
			}
		}
	}

	if (!has_acid and !has_basic and !has_polar and !has_all){
		for (unsigned j=0; j<resnames.size(); j++){
			this->mol_extraction(M1, MExtract1, resnames[j]);
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
            if ((resnames[i] == unique[j])){
				found = true;
			}
		}
		if (!found){
			unique.push_back(resnames[i]);
		}
	}
	return(unique);
}

int Engine::serial_search(Mol* ME1, Parser* Input, vector<string> unique, vector<string> pdb_list) {
	for (unsigned i=0; i< pdb_list.size(); i++) {
        Mol* M2 = new Mol(Input);
        Mol* MExtract2 = new Mol(Input);
		if (M2->read_pdb(pdb_list[i])){
			this->mol_extraction(M2, MExtract2, unique);

			if (Input->verbose){
                sprintf(info, "Found %d matching residues in search molecule.", int(MExtract2->mymol.size()));
                Writer->print_info(info);
				for (unsigned a=0; a<MExtract2->mymol.size(); a++){
                    sprintf(info, "        %s%d", MExtract2->mymol[a].resname.c_str(), MExtract2->mymol[a].resnumber);
                    Writer->print_info(info);
				}
			}


			vector<vector<vector<double> > >xyz;
            Optimization* Opt = new Optimization(Writer, ME1, Input);
			opt_result_t* opt_result = new opt_result_t;
            opt_result->succeded = false;
            if (MExtract2->mymol.size() > 0){
                Opt->optimize_rmsd(MExtract2, opt_result);
            }

            if (opt_result->succeded){
                sprintf(info, "    Matched residues:");
                Writer->print_info(info);
                for (unsigned a=0; a< opt_result->imatched1.size(); a++){
                    sprintf(info, "        %15.15s %4.4s %4d       %15.15s %4.4s %4d %6.4f", ME1->filename.c_str(), opt_result->smatched1[a].c_str(), opt_result->imatched1[a],
                           MExtract2->filename.c_str(), opt_result->smatched2[a].c_str(), opt_result->imatched2[a], opt_result->rmsds[a]);
                    Writer->print_info(info);
                }
            }

            else {
                if (Input->verbose){
                    printf("No matched results for %-15.15s", MExtract2->filename.c_str());
                }
            }

            if (opt_result->succeded && Input->write_pdb){
                Coord* CoordManip = new Coord;
                xyz = CoordManip->rototranslate(M2, opt_result->rotrans[0], opt_result->rotrans[1], opt_result->rotrans[2], opt_result->rotrans[3],
                        opt_result->rotrans[4], opt_result->rotrans[5]);
                delete CoordManip;
                if (Input->verbose){
                    sprintf(info, "Writting pdb %s ....", (pdb_list[i].substr(0,pdb_list[i].find(".pdb")) + "_smatch").c_str());
                    Writer->print_info(info);
                }
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
int Engine::run_over_mpi(int argc, char* argv[], Mol* ME1, vector<string> unique, Parser* Input){
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

        if (Input->verbose){
            sprintf(info, "There are %d files, divided into %d chuncks of size %d.", int(pdb_list.size()), world.size(), chunck_size);
            Writer->print_info(info);
        }

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

        if (Input->verbose){
            for (unsigned i=0; i< chuncks.size(); i++){
                sprintf(info, "Chuncks[%d] has %d elements.", int(i), int(chuncks[i].size()));
                Writer->print_info(info);
            }
        }

		scatter(world, chuncks, tmp,0);
        serial_search(ME1, Input, unique, tmp);
	}
	else{
		scatter(world, tmp, 0);
        serial_search(ME1, Input, unique, tmp);
	}
	return 0;
}
#endif

int Engine::run_serial(Mol* ME1, vector<string> unique, Parser* Input){
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

    this->serial_search_omp(ME1, Input, unique, pdb_list);

	return 0;

}

int Engine::serial_search_omp(Mol* ME1, Parser* Input, vector<string> unique, vector<string> pdb_list) {

#pragma omp parallel num_threads(Input->parallel_jobs)

	{
#pragma omp for schedule(static, 1)
		for (unsigned i=0; i< pdb_list.size(); i++) {
            Mol* M2 = new Mol(Input);
            Mol* MExtract2 = new Mol(Input);
			if (M2->read_pdb(pdb_list[i])){
				this->mol_extraction(M2, MExtract2, unique);

                if (Input->verbose){
                    sprintf(info, "Found %d matching residues in search molecule.", int(MExtract2->mymol.size()));
                    Writer->print_info(info);
                }

				vector<vector<vector<double> > >xyz;
                opt_result_t* opt_result = new opt_result_t;
                opt_result->succeded = false;
                Optimization* Opt = new Optimization(Writer, ME1, Input);
                Opt->optimize_rmsd(MExtract2, opt_result);

                if (opt_result->succeded){
                    sprintf(info, "    Matched residues:");
                    Writer->print_info(info);
                    for (unsigned a=0; a< opt_result->imatched1.size(); a++){
                        sprintf(info, "        %15.15s %4.4s %4d       %15.15s %4.4s %4d %6.4f", ME1->filename.c_str(), opt_result->smatched1[a].c_str(), opt_result->imatched1[a],
                               MExtract2->filename.c_str(), opt_result->smatched2[a].c_str(), opt_result->imatched2[a], opt_result->rmsds[a]);
                        Writer->print_info(info);
                    }
                }

                else {
                    if (Input->verbose){
                        sprintf(info, "No matched results for %-15.15s", MExtract2->filename.c_str());
                        Writer->print_info(info);
                    }
                }

				if (opt_result->succeded and Input->write_pdb){

                    Coord* CoordManip = new Coord;
                    xyz.clear();
                    xyz = CoordManip->rototranslate_new_ref(M2, MExtract2, opt_result->rotrans[0], opt_result->rotrans[1], opt_result->rotrans[2], opt_result->rotrans[3],
							opt_result->rotrans[4], opt_result->rotrans[5]);
					delete CoordManip;
					Writer->write_pdb(M2, xyz, 0.0, opt_result->rmsd, (pdb_list[i].substr(0,pdb_list[i].find(".pdb")) + "_smatch"));
//                    Writer->write_pdb(MExtract2, opt_result->xyz, 0.0, opt_result->rmsd, (pdb_list[i].substr(0,pdb_list[i].find(".pdb")) + "_smatch_site"));
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

bool Engine::residues_match(string r1, string r2){
	bool result = false;

	if (r1 == r2){
		result = true;
	}

	else if (r2 == "acid" and ((r1 == "ASP") or (r1 == "GLU"))){
		result = true;
	}

	else if (r2 == "basic" and ((r1 == "LYS") or (r1 == "ARG"))){
		result = true;
	}

	else if (r2 == "*"){
		result = true;
	}
	return result;
}
