/*
 * SMatch.h
 *
 *  Created on: Oct 22, 2013
 *      Author: asn
 */

#ifndef SMATCH_H_
#define SMATCH_H_

struct opt_result_t {
		bool succeded;
		vector<double> rotrans;
	    vector<vector<vector<double> > >xyz;
		double rmsd;
	};

#endif /* SMATCH_H_ */
