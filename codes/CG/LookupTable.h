/*
 * LookupTable.h
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#ifndef CODES_LOOKUPTABLE_H_
#define CODES_LOOKUPTABLE_H_

#include "definitions.h"

#include "../utils/strings.h"

#include <vector>
#include <string>
#include <iostream>

class LookupTable;

class Mesh {
public:
	Mesh() {
		delta = inv_sqr_delta = xlow = xupp = -1.;
	};

	void build(LookupTable *that, double (LookupTable::*f)(double, void*), double (LookupTable::*der)(double, void*), void *args, int npoints, double mxlow, double mxupp);

	double query_derivative(double x);

	double query(double x);

	~Mesh() {

	}

	double delta, inv_sqr_delta, xlow, xupp;
	std::vector<double> A, B, C, D;
};

class LookupTable {
	struct base_function {
		std::vector<double> x, fx, dfx;

		bool add_line(std::string line) {
			std::vector<std::string> spl = split(line, ' ');
			if(spl.size() != 3) {
				return false;
			}

			double new_x = atof(spl[0].c_str());
			double new_fx = atof(spl[1].c_str());
			double new_dfx = atof(spl[2].c_str());

			x.push_back(new_x);
			fx.push_back(new_fx);
			dfx.push_back(new_dfx);

			return true;
		}
	};

public:
	LookupTable(std::string filename, int points);
	virtual ~LookupTable();

	void potential(CGBead &p, CGBead &q, double shift_by);

private:
	double _linear_interpolation(double x, std::vector<double> &x_data, std::vector<double> &fx_data);

	double _fx(double x, void *par);

	double _dfx(double x, void *par);

	Mesh _lookup_table;
	int _lt_points;
};

#endif /* CODES_LOOKUPTABLE_H_ */
