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
	return (xyz);
}

double Optimization::pre_optimize_rmsd_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
	opt_data* odata = (opt_data*) data;
	Coord* Manip = new Coord;
    vector<vector<vector<double> > >new_xyz;
	double f;

    new_xyz = Manip->rototranslate(odata->M2, x[0], x[1], x[2], x[3], x[4], x[5]);
	f = Optimization::compute_rmsd(M1, odata->M2, new_xyz, 0, odata->resnumber);
	return(f);
}

void Optimization::optimize_rmsd(Mol* M2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LN_NELDERMEAD,6);
	opt_data *odata = new opt_data;
    vector<vector<vector<double> > >xyz;
	odata->M2 = M2;
	opt_result->succeded = false;

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
	double rmsd_total, rmsd, optimal_rmsd;
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
//			printf("Comparing residue: M2[%d] (%.3f)\n", i, fo);
            for (unsigned k=0; k< M1->mymol.size(); k++){
                for (unsigned j=0; j<M2->mymol.size(); j++){
                    if (M1->mymol[k].resname == M2->mymol[j].resname){
						rmsd = this->compute_rmsd(M1, M2, xyz, k, j);
//						printf("Alignment of M1[%d] with M2[%d] = %.3f\n", k, j, rmsd);
						if (rmsd <= 5.0){
							rmsd_total+= rmsd;
							nres_sol++;
						}
					}
				}
			}

//			printf("RMSD: %.3f for %d residues\n", rmsd_total, nres_sol);

            if ((rmsd_total < optimal_rmsd) and (nres_sol == int(M1->mymol.size()))){
				optimal_rmsd=rmsd_total;
				optimal_x = x;
				optimal_Nres=nres_sol;
			}
		}
	}
#ifdef DEBUG
	printf("alpha: %8.3f beta: %8.3f gamma: %8.3f\n", optimal_x[0], optimal_x[1], optimal_x[2]);
	printf("x: %8.3f y:%8.3f z:%8.3f\n", optimal_x[3], optimal_x[4], optimal_x[5]);
#endif
	xyz = update_coords(optimal_x, M2);
	if (optimal_rmsd < 50.0){
		sprintf(info, "FILE = %-40.40s RMSD = %8.3f  N = %4d",M2->filename.c_str(), optimal_rmsd, optimal_Nres);
		Writer->print_info(info);

		opt_result->rmsd = optimal_rmsd;
		opt_result->xyz = xyz;
		opt_result->rotrans = optimal_x;
		opt_result->succeded = true;
	}
	delete opt;
}
