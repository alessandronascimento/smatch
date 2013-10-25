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
        		printf("ATOM   %4d %-3s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f  0.00    0.00  1\n", atom_count, Cmol->mymol[i].atomnames[j].c_str(),
        				Cmol->mymol[i].resname.c_str(),  j+1, Cmol->mymol[i].xyz[j][0], Cmol->mymol[i].xyz[j][1], Cmol->mymol[i].xyz[j][2]);
        		gzprintf(outpdb, "ATOM   %4d %-3s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f  0.00    0.00  1\n", atom_count, Cmol->mymol[i].atomnames[j].c_str(),
        				Cmol->mymol[i].resname.c_str(),  j+1, Cmol->mymol[i].xyz[j][0], Cmol->mymol[i].xyz[j][1], Cmol->mymol[i].xyz[j][2]);
        		atom_count++;
        	}
        }
        gzprintf(outpdb, "TER\n");
        gzclose(outpdb);
}

void Printer::print_info(char info[98]){
	fprintf(logfile, "%-98s\n", info);
}
