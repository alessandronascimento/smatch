/*
 * Gaussian.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: asn
 */

#include "Gaussian.h"

double Gaussian::dist_squared(double x1, double x2, double y1, double y2, double z1, double z2){
	double r2 = ((x2-x1)*(x2-x1)) + ((y2-y1)*(y2-y1)) + ((z2-z1)*(z2-z1));
	return (r2);
}

Gaussian::Gaussian() {
	// TODO Auto-generated constructor stub
}

Gaussian::~Gaussian() {
	// TODO Auto-generated destructor stub
}

double Gaussian::shape_and_charge_density(Mol* M1, Mol* M2){
	double pi=2.0*sqrt(2);
	double pj=2.0*sqrt(2);
	double dij2;
	double Vshape=0.0, Vpos=0.0, Vneg=0.0;
	double alphai_shape=0.0, alphai_pos=0.0, alphai_neg=0.0;
	double alphaj_shape=0.0, alphaj_pos=0.0, alphaj_neg=0.0;
	for (int i=0; i<M1->N; i++){
		if (M1->radii[i] != 0.0){
			alphai_shape = PI*( pow( ( (3*pi)/(4*PI*pow(M1->radii[i], 3))), (2.0/3.0)));
		}
		else {
			alphai_shape=0.0;
		}

		if (M1->charges[i] > 0.0){
			alphai_pos = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + M1->charges[i]), 3))), (2.0/3.0)));
			alphai_neg=0.0;
		}

		else if (M1->charges[i] < 0.0) {
			alphai_neg = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + abs(M1->charges[i])), 3))), (2.0/3.0)));
			alphai_pos=0.0;
		}

		for (int j=0; j<M2->N; j++){
			if (M2->radii[j] != 0.0000){
				alphaj_shape = PI*pow(  ( (3*pi)/(4*PI*pow(M2->radii[j], 3))), (2.0/3.0));
			}
			else {
				alphaj_shape=0.0000;
			}

			if (M2->charges[j] > 0.00){
				alphaj_pos = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + M2->charges[j]), 3))), (2.0/3.0));
				alphaj_neg=0.0;
			}

			else if (M2->charges[j] < 0.00){
				alphaj_neg = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + abs(M2->charges[j])), 3))), (2.0/3.0));
				alphaj_pos = 0.0;
			}

			dij2 = this->dist_squared(M1->xyz[i][0], M2->xyz[j][0],M1->xyz[i][1], M2->xyz[j][1],M1->xyz[i][2], M2->xyz[j][2]);

			if (alphai_shape != 0.0 and alphaj_shape != 0.0 ){
				Vshape += pi*pj*(exp((-alphai_shape*alphaj_shape*dij2)/(alphai_shape+alphaj_shape)))*(pow((PI/(alphai_shape+alphaj_shape)), (3.0/2.0)));
			}

			if (alphai_pos != 0.0 and alphaj_pos != 0.0 ){
				Vpos += pi*pj*(exp((-alphai_pos*alphaj_pos*dij2)/(alphai_pos+alphaj_pos)))*(pow((PI/(alphai_pos+alphaj_pos)), (3.0/2.0)));
			}

			if (alphai_neg != 0.0 and alphaj_neg != 0.0 ){
				Vneg += pi*pj*(exp((-alphai_neg*alphaj_neg*dij2)/(alphai_neg+alphaj_neg)))*(pow((PI/(alphai_neg+alphaj_neg)), (3.0/2.0)));
			}
		}
	}
	return(Vshape+Vpos+Vneg);
}
