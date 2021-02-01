/*
 * LookupTable.h
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#ifndef CODES_LOOKUPTABLE_H_
#define CODES_LOOKUPTABLE_H_

#include "definitions.h"

#include <vector>
#include <string>

class LookupTable;

class Mesh {
public:
	Mesh() : N(0) {
		delta = inv_sqr_delta = xlow = xupp = -1.;
	};

	void init(int size) {
		N = size;
		delta = 0;
		A.resize(size + 1);
		B.resize(size + 1);
		C.resize(size + 1);
		D.resize(size + 1);
	}

	void build(LookupTable *that, double (LookupTable::*f)(double, void*), double (LookupTable::*der)(double, void*), void *args, int npoints, double mxlow, double mxupp);

	double query_derivative(double x);

	double query(double x);

	~Mesh() {

	}

	int N;
	double delta, inv_sqr_delta, xlow, xupp;
	std::vector<double> A, B, C, D;
};

class LookupTable {
	struct base_function {
		int points;
		std::vector<double> x, fx, dfx;
	};

public:
	LookupTable(std::string filename, int points);
	virtual ~LookupTable();

	void potential(Bead &p, Bead &q, double shift_by);

private:
	double _linear_interpolation(double x, std::vector<double> &x_data, std::vector<double> fx_data, int points);

	double _fx(double x, void *par);

	double _dfx(double x, void *par);

	Mesh _lookup_table;
	int _lt_points;
};

#endif /* CODES_LOOKUPTABLE_H_ */
