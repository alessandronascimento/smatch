/*
 * Printer.cpp
 *
 *  Created on: Oct 2, 2013
 *      Author: asn
 */

#include "Printer.h"

Printer::Printer() {
	logfile = fopen("SMatch.log", "w");
}

Printer::Printer(Parser* Input){
	logfile = fopen((Input->output_prefix + ".log").c_str(), "w");
    this->print_params(Input);
}

Printer::~Printer() {
	fclose(logfile);
}

void Printer::write_pdb(Mol *Cmol, double energy, double rmsd, string outname){
        gzFile outpdb;
        outpdb = gzopen((outname+".pdb.gz").c_str(), "w");
      gzprintf(outpdb, "REMARK\n");
      gzprintf(outpdb, "REMARK %-9s energy = %9.2f rmsd   = %12.3f\n", outname.c_str(), energy, rmsd);
        int atom_count=1;
        for (unsigned i=0; i<Cmol->mymol.size(); i++){
        	for (unsigned j=0; j< Cmol->mymol[i].atomnames.size(); j++){
        		gzprintf(outpdb, "ATOM   %4d %-3s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f  0.00    0.00  1\n", atom_count, Cmol->mymol[i].atomnames[j].c_str(),
                        Cmol->mymol[i].resname.c_str(),  Cmol->mymol[i].resnumber, Cmol->mymol[i].xyz[j][0], Cmol->mymol[i].xyz[j][1], Cmol->mymol[i].xyz[j][2]);
        		atom_count++;
        	}
        }
        gzprintf(outpdb, "TER\n");
        gzclose(outpdb);
}

void Printer::write_pdb(Mol *Cmol, vector<vector<vector<double> > > xyz, double energy, double rmsd, string outname){
        gzFile outpdb;
        outpdb = gzopen((outname+".pdb.gz").c_str(), "w");
      gzprintf(outpdb, "REMARK\n");
      gzprintf(outpdb, "REMARK %-9s energy = %9.2f rmsd   = %12.3f\n", outname.c_str(), energy, rmsd);
        int atom_count=1;
        for (unsigned i=0; i<Cmol->mymol.size(); i++){
        	for (unsigned j=0; j< Cmol->mymol[i].atomnames.size(); j++){
        		gzprintf(outpdb, "ATOM   %4d %-3s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f  0.00    0.00  1\n", atom_count, Cmol->mymol[i].atomnames[j].c_str(),
                        Cmol->mymol[i].resname.c_str(),  Cmol->mymol[i].resnumber, xyz[i][j][0], xyz[i][j][1], xyz[i][j][2]);
        		atom_count++;
        	}
        }
        gzprintf(outpdb, "TER\n");
        gzclose(outpdb);
}

void Printer::print_info(char info[98]){
	fprintf(logfile, "%-98s\n", info);
	printf("%s\n", info);
}

void Printer::print_welcome(void){
    printf("****************************************************************************\n");
    printf("*            SMATCH - Structure-Based Active Site Matching                 *\n");
    printf("*                                                                          *\n");
    printf("*  Alessandro S. Nascimento - Instituto de Fisica de Sao Carlos - IFSC/USP *\n");
    printf("*                         asnascimento@ifsc.usp.br                         *\n");
    printf("*                    http://www.biotechmol.ifsc.usp.br                     *\n");
    printf("*                                                                          *\n");
    printf("****************************************************************************\n");
}

void Printer::print_params(Parser* Input){
    printf("****************************************************************************\n");
    printf("* %20.20s               %30.30s*\n", "reference_file", Input->reference_file.c_str());
    printf("* %20.20s               %30f*\n", "search radius", Input->search_radius);
    printf("* %20.20s               %30.30s*\n", "output_prefix", Input->output_prefix.c_str());
    printf("* %20.20s               %30d*\n", "write_pdb", Input->write_pdb);
    printf("* %20.20s               %30d*\n", "matching_residues", Input->matching_residues);
    printf("* %20.20s               %30.30s*\n", "directory", Input->directory.c_str());
    printf("****************************************************************************\n");
}
