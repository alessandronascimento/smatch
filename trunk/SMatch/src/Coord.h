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
	vector<double> compute_com(vector<vector<double> >coords, Mol *M1);
	vector<vector<double> > translate(Mol* M1, double dx, double dy, double dz);
	vector<vector<double> >rototranslate(vector<vector<double> >coordinates, Mol* M1, double alpha, double beta, double gamma, double transx, double transy, double transz);
	virtual ~Coord();
};

#endif /* COORD_H_ */
