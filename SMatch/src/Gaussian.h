/*
 * Gaussian.h
 *
 *  Created on: Oct 1, 2013
 *      Author: asn
 */

#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include "Mol.h"
#include <cmath>

#define PI 3.14159265359

using namespace std;

class Gaussian {
public:
	Gaussian();
	double shape_and_charge_density(Mol* M1, Mol* M2);
	double dist_squared(double x1, double x2, double y1, double y2, double z1, double z2);
	virtual ~Gaussian();
};

#endif /* GAUSSIAN_H_ */
