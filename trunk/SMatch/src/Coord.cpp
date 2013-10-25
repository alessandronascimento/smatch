/*
 * Coord.cpp
 *
 *  Created on: Oct 2, 2013
 *      Author: asn
 */

#include "Coord.h"

Coord::Coord() {
	// TODO Auto-generated constructor stub

}

Coord::~Coord() {
	// TODO Auto-generated destructor stub
}


vector<double> Coord::compute_com(vector<vector<double> >coords, Mol *M1){
	double centerx=0.0;
	double centery=0.0;
	double centerz=0.0;
	double totalmass=0.0;
	vector<double> com(3);
	for(int i=0; i<M1->N; i++){
		centerx+= (12.0*coords[i][0]);
		centery+= (12.0*coords[i][1]);
		centerz+= (12.0*coords[i][2]);
		totalmass+= 12.0;
	}
	com[0] = (centerx/totalmass);
	com[1] = (centery/totalmass);
	com[2] = (centerz/totalmass);
	return(com);
}

vector<vector<double> >Coord::rototranslate(vector<vector<double> >coordinates, Mol* M1, double alpha, double beta, double gamma, double transx, double transy, double transz){
	vector<vector<double> >new_coordinates;
	vector<double> txyz(3);
	vector<double> COM = this->compute_com(coordinates, M1);
	double x, y, z;
	for(int i=0; i < M1->N ; i++){
		x=coordinates[i][0]-COM[0];
		y=coordinates[i][1]-COM[1];
		z=coordinates[i][2]-COM[2];
		txyz[0] = ((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx+COM[0]);
		txyz[1] = ((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy + COM[1]);
		txyz[2] = ((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz + COM[2]);
		new_coordinates.push_back(txyz);
	}
	return(new_coordinates);
}

vector<vector<double> > Coord::translate(Mol* M1, double dx, double dy, double dz){
	vector<vector<double> > new_coordinates;
	vector<double> t;
	for (int i=0; i< M1->N; i++){
		t.push_back(M1->xyz[i][0]+dx);
		t.push_back(M1->xyz[i][1]+dy);
		t.push_back(M1->xyz[i][2]+dz);
		new_coordinates.push_back(t);
		t.clear();
	}
	return(new_coordinates);
}
