/*
 * Printer.h
 *
 *  Created on: Oct 2, 2013
 *      Author: asn
 */

#ifndef PRINTER_H_
#define PRINTER_H_

#include <zlib.h>
#include <string>
#include "Mol.h"
#include "Parser.h"

using namespace std;

class Printer {
public:
	FILE* logfile;
	Printer();
	Printer(Parser* Input);
	void write_pdb(Mol *Cmol, double energy, double rmsd, string outname);
	void write_pdb(Mol *Cmol, vector<vector<vector<double> > > xyz, double energy, double rmsd, string outname);
	virtual ~Printer();
	void print_info(char info[98]);
};

#endif /* PRINTER_H_ */
