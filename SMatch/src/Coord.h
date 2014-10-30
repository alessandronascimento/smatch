/*
 * Coord.h
 *
 *  Created on: Oct 2, 2013
 *      Author: asn
 */

#ifndef COORD_H_
#define COORD_H_

#include <vector>
#include <string>
#include <cmath>
#include "Mol.h"

#define PI 3.14159265359

using namespace std;

class Coord {
public:
	Coord();
    vector<double> compute_com(Mol *M1);
    vector<vector<vector<double> > >rototranslate(Mol* M1, double alpha, double beta, double gamma, double transx, double transy, double transz);
    vector<vector<vector<double> > >rototranslate_new_ref(Mol* M1, Mol* M2, double alpha, double beta, double gamma, double transx, double transy, double transz);
	virtual ~Coord();
};

#endif /* COORD_H_ */
