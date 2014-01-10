/*
 * Optimization.cpp
 *
 *  Created on: Oct 2, 2013
 *      Author: asn
 */

#include "Optimization.h"

Optimization::Optimization(Printer* _Writer, Mol* _M1) {
	M1 = _M1;
	Writer = _Writer;

}

Optimization::Optimization(Printer* _Writer, Mol* _M1, Parser* _Input) {
    M1 = _M1;
    Writer = _Writer;
    Input = _Input;
}

Optimization::~Optimization() {
}

double Optimization::dist_squared(double x1, double x2, double y1, double y2, double z1, double z2){
	return(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1)));
}

double Optimization::compute_rmsd(Mol* M1, Mol* M2, vector<vector<vector<double> > >xyz, int r1, int r2){
	double rmsd=0.0;
	int Natoms=0;
    for (unsigned i=0; i<M1->mymol[r1].atomnames.size(); i++){
        for (unsigned j=0; j<M2->mymol[r2].atomnames.size(); j++){
            if (M1->mymol[r1].atomnames[i] == M2->mymol[r2].atomnames[j]){
                rmsd += dist_squared(M1->mymol[r1].xyz[i][0], xyz[r2][j][0], M1->mymol[r1].xyz[i][1],
                        xyz[r2][j][1],M1->mymol[r1].xyz[i][2], xyz[r2][j][2]);
                Natoms++;
            }
        }
    }
	rmsd=(rmsd/Natoms);
	rmsd=sqrt(rmsd);
	return(rmsd);
}

vector<vector<vector<double> > >Optimization::update_coords(const std::vector<double> &x, Mol* M2){
    vector<vector<vector<double> > >xyz;
	Coord* Manip = new Coord;
	if (x.size() == 6){
        xyz = Manip->rototranslate(M2, x[0], x[1], x[2], x[3], x[4], x[5]);
	}
	else{
        for (unsigned i=0; i< M2->mymol.size(); i++){       // copying M2->mymol->xyz to
            xyz.push_back(M2->mymol[i].xyz);                // xyz.
        }
	}
	delete Manip;
	return (xyz);
}

double Optimization::pre_optimize_rmsd_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
	opt_data* odata = (opt_data*) data;
	Coord* Manip = new Coord;
    vector<vector<vector<double> > >new_xyz;
	double f;

    new_xyz = Manip->rototranslate(odata->M2, x[0], x[1], x[2], x[3], x[4], x[5]);
	f = Optimization::compute_rmsd(M1, odata->M2, new_xyz, 0, odata->resnumber);
	delete Manip;
	return(f);
}

void Optimization::optimize_rmsd(Mol* M2, opt_result_t* opt_result){
    int nmatch = Input->matching_residues;
    if ((nmatch < 1) or (nmatch > int(M1->mymol.size()))){
        nmatch = int(M1->mymol.size());
    }

    nlopt::opt *opt = new nlopt::opt(nlopt::LN_NELDERMEAD,6);
	opt_data *odata = new opt_data;
    vector<vector<vector<double> > >xyz;
	odata->M2 = M2;
	opt_result->succeded = false;

    vector<int> imatched1, imatched2;
    vector<string> smatched1, smatched2;
    vector<double> rmsds;

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = -50.0;
	lb[4] = -50.0;
	lb[5] = -50.0;
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = 50.0;
	ub[4] = 50.0;
	ub[5] = 50.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);
	opt->set_xtol_rel(1.0E-10);
	opt->set_maxtime(60);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	opt->set_min_objective(Optimization::pre_optimize_rmsd_function, odata);

	double fo;
    double rmsd_total, rmsd=-1.0, optimal_rmsd;
    vector<double> optimal_x=x;
    int optimal_Nres=0;
    optimal_rmsd=9999.0;

    for (unsigned i=0; i<M2->mymol.size(); i++){
    	if (M2->mymol[i].resname == M1->mymol[0].resname){
    		x[0] = 0.0;
    		x[1] = 0.0;
    		x[2] = 0.0;
    		x[3] = 0.0;
    		x[4] = 0.0;
    		x[5] = 0.0;
    		odata->resnumber = int(i);
    		opt->optimize(x,fo);
    		xyz = update_coords(x, M2);
    		rmsd_total = 0.00;
    		int nres_sol=0;
    		for (unsigned k=0; k< M1->mymol.size(); k++){
    			for (unsigned j=0; j<M2->mymol.size(); j++){
    				if (this->residues_match(M2->mymol[j].resname, Input->lookup[k])){
    					if ((M1->mymol[k].resname == M2->mymol[j].resname)){
    						rmsd = this->compute_rmsd(M1, M2, xyz, k, j);
    					}
    					else {
    						rmsd = this->compute_rmsd_non_similar(M1, M2, k, j);
    					}
    					if ((rmsd >= 0.0) and (rmsd <= Input->search_radius)){
    						rmsd_total+= rmsd;
    						nres_sol++;
    						smatched1.push_back(M1->mymol[k].resname);
    						smatched2.push_back(M2->mymol[j].resname);
    						imatched1.push_back(M1->mymol[k].resnumber);
    						imatched2.push_back(M2->mymol[j].resnumber);
    						rmsds.push_back(rmsd);
    					}
    				}
    			}
    		}
    		if ((rmsd_total < optimal_rmsd) and (nres_sol >= nmatch)){
    			optimal_rmsd=rmsd_total;
    			optimal_x = x;
    			optimal_Nres=nres_sol;
    			opt_result->imatched1 = imatched1;
    			opt_result->imatched2 = imatched2;
    			opt_result->smatched1 = smatched1;
    			opt_result->smatched2 = smatched2;
    			opt_result->rmsds = rmsds;
    		}
    	}
    	imatched1.clear();
    	imatched2.clear();
    	smatched1.clear();
    	smatched2.clear();
    	rmsds.clear();
    }

    xyz = update_coords(optimal_x, M2);
    if (optimal_rmsd < (2.0*nmatch*Input->search_radius)){
    	sprintf(info, "FILE = %-40.40s RMSD = %8.3f  N = %4d",M2->filename.c_str(), optimal_rmsd, optimal_Nres);
    	Writer->print_info(info);

    	opt_result->rmsd = optimal_rmsd;
    	opt_result->xyz = xyz;
    	opt_result->rotrans = optimal_x;
    	opt_result->succeded = true;
    }
    delete opt;
    delete odata;
}

double Optimization::compute_rmsd_non_similar(Mol* M1, Mol* M2, int r1, int r2){
	double rmsd=0.0;
	int natom=0;
	for (unsigned i=0; i< M1->mymol[r1].xyz.size(); i++){
		for (unsigned j=0; j < M2->mymol[r2].xyz.size(); j++){
			if (M1->mymol[r1].atomnames[i] == M2->mymol[r2].atomnames[j]){
				rmsd += this->dist_squared(M1->mymol[r1].xyz[i][0], M2->mymol[r2].xyz[j][0], M1->mymol[r1].xyz[i][1], M2->mymol[r2].xyz[j][1], M1->mymol[r1].xyz[i][2], M2->mymol[r2].xyz[j][2]);
				natom++;
			}
		}
	}
	if (natom != 0){
		rmsd = rmsd/natom;
		rmsd = sqrt(rmsd);
	}
	else {
		rmsd = -1.0;
	}
	return (rmsd);
}

bool Optimization::residues_match(string r1, string r2){
	bool result = false;

	if (r1 == r2){
		result = true;
	}

	else if (r2 == "acid" and ((r1 == "ASP") or (r1 == "GLU"))){
		result = true;
	}

	else if (r2 == "basic" and ((r1 == "LYS") or (r1 == "ARG"))){
		result = true;
	}

	else if (r2 == "*"){
		result = true;
	}
	return result;
}
