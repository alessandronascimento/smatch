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

vector<double> Coord::compute_com(Mol *M1){
    double centerx=0.0;
    double centery=0.0;
    double centerz=0.0;
    double totalmass=0.0;
    vector<double> com(3);
    for(unsigned i=0; i<M1->mymol.size(); i++){
        for (unsigned j=0; j<M1->mymol[i].xyz.size(); j++){
            centerx+= (12.0*M1->mymol[i].xyz[j][0]);
            centery+= (12.0*M1->mymol[i].xyz[j][1]);
            centerz+= (12.0*M1->mymol[i].xyz[j][2]);
            totalmass+= 12.0;
        }
    }
    com[0] = (centerx/totalmass);
    com[1] = (centery/totalmass);
    com[2] = (centerz/totalmass);
    return(com);
}

vector<vector<vector<double> > >Coord::rototranslate(Mol* M1, double alpha, double beta, double gamma, double transx, double transy, double transz){
    vector<vector<double> >nc;
    vector<vector<vector<double> > >new_coordinates;
    vector<double> txyz(3);
    vector<double> COM = this->compute_com(M1);
    double x, y, z;
    for (unsigned i=0; i<M1->mymol.size(); i++){
        for (unsigned j=0; j<M1->mymol[i].xyz.size(); j++){
            x = M1->mymol[i].xyz[j][0] - COM[0];
            y = M1->mymol[i].xyz[j][1] - COM[1];
            z = M1->mymol[i].xyz[j][2] - COM[2];

            txyz[0] = ((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+ transx + COM[0]);
            txyz[1] = ((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+ transy + COM[1]);
            txyz[2] = ((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+ transz + COM[2]);
            nc.push_back(txyz);
        }
        new_coordinates.push_back(nc);
        nc.clear();
    }
    return(new_coordinates);
}


vector<vector<vector<double> > >Coord::rototranslate_new_ref(Mol* M1, Mol* M2, double alpha, double beta, double gamma, double transx, double transy, double transz){
    vector<vector<double> >nc;
    vector<vector<vector<double> > >new_coordinates;
    vector<double> txyz(3);
    vector<double> COM = this->compute_com(M2);
    double x, y, z;
    for (unsigned i=0; i<M1->mymol.size(); i++){
        for (unsigned j=0; j<M1->mymol[i].xyz.size(); j++){
            x = M1->mymol[i].xyz[j][0] - COM[0];
            y = M1->mymol[i].xyz[j][1] - COM[1];
            z = M1->mymol[i].xyz[j][2] - COM[2];

            txyz[0] = ((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+ transx + COM[0]);
            txyz[1] = ((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+ transy + COM[1]);
            txyz[2] = ((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+ transz + COM[2]);
            nc.push_back(txyz);
        }
        new_coordinates.push_back(nc);
        nc.clear();
    }
    return(new_coordinates);
}
