/*
 * Optimization.h
 *
 *  Created on: Oct 2, 2013
 *      Author: asn
 */

#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

#include <vector>
#include <nlopt.h>
#include <nlopt.hpp>
#include "Mol.h"
#include "Parser.h"
#include "Coord.h"
#include "Printer.h"
#include "SMatch.h"
#include "Parser.h"

class Optimization {
public:

	static Mol* M1;
	Printer* Writer;
    Parser* Input;
	char info[98];
	Optimization(Printer* _Writer, Mol* _M1);
    Optimization(Printer* _Writer, Mol* _M1, Parser* _Input);
	static double dist_squared(double x1, double x2, double y1, double y2, double z1, double z2);
	static double pre_optimize_rmsd_function(const std::vector<double> &x, std::vector<double> &grad, void *data);
	static double post_optimize_rmsd_function(const std::vector<double> &x, std::vector<double> &grad, void *data);
	void minimize_overlay_nlopt_ln_auglag(Mol* M2);
    vector<vector<vector<double> > >update_coords(const std::vector<double> &x, Mol* M2);
	virtual ~Optimization();
    static double compute_rmsd(Mol* M1, Mol* M2, vector<vector<vector<double> > >xyz, int r1, int r2);
	void optimize_rmsd(Mol* M2, opt_result_t* opt_result);
    double compute_rmsd_non_similar(Mol* M1, Mol* M2, int r1, int r2);
    static double compute_rmsd_non_similar(Mol* M1, Mol* M2, vector<vector<vector<double> > > xyz, int r1, int r2);
    bool residues_match(string r1, string r2);
	struct opt_data{
		Mol* M2;
		int resnumber;
	};

	struct opt_vector_data{
		Mol* M2;
		vector<int> m1_residues;
		vector<int> m2_residues;
	};
};

#endif /* OPTIMIZATION_H_ */
