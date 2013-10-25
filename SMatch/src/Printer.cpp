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

void Printer::write_pdb(Mol *Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname){
        gzFile outpdb;
        outpdb = gzopen((outname+".pdb.gz").c_str(), "w");
//      gzprintf(outpdb, "MDL\n");
      gzprintf(outpdb, "REMARK\n");
      gzprintf(outpdb, "REMARK %-9s energy = %9.2f rmsd   = %12.3f\n", outname.c_str(), energy, rmsd);
        int i=0;
        int resn=0;

        while (resn < int(Cmol->residue_pointer.size()-1)){
                while(i < Cmol->residue_pointer[resn+1]){
                        gzprintf(outpdb, "ATOM   %4d %-3s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f  0.00    0.00  1\n", i, Cmol->atomnames[i].c_str(), Cmol->resnames[resn].c_str(), resn+1, xyz[i][0], xyz[i][1], xyz[i][2]);
                        i++;
                }
                resn++;
        }
        while (i < Cmol->N){
                gzprintf(outpdb, "ATOM   %4d %-3s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f  0.00    0.00  1\n", i, Cmol->atomnames[i].c_str(), Cmol->resnames[resn].c_str(), resn+1, xyz[i][0], xyz[i][1], xyz[i][2]);
                i++;
        }
        gzprintf(outpdb, "TER\n");
        gzprintf(outpdb, "ENDMDL\n");
        gzclose(outpdb);
}

void Printer::print_info(char info[98]){
	fprintf(logfile, "%-98s\n", info);
}
