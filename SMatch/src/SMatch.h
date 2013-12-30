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
        vector<string> smatched1;
        vector<string> smatched2;
        vector<int> imatched1;
        vector<int> imatched2;
        vector<double> rmsds;
	};

#endif /* SMATCH_H_ */
