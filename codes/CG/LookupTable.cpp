/*
 * LookupTable.cpp
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#include "LookupTable.h"

#include "../utils/strings.h"

#include <iostream>
#include <fstream>
#include <cfloat>

void Mesh::build(LookupTable *that, double (LookupTable::*f)(double, void*), double (LookupTable::*der)(double, void*), void *args, int npoints, double mxlow, double mxupp) {
	int i;
	double x;

	init(npoints);

	double dx = (xupp - xlow) / (double) npoints;
	delta = dx;
	inv_sqr_delta = 1 / (dx * dx);
	xlow = mxlow;
	xupp = mxupp;

	double fx0, fx1, derx0, derx1;

	for(i = 0; i < npoints + 1; i++) {
		x = xlow + i * dx;

		fx0 = (that->*f)(x, args);
		fx1 = (that->*f)(x + dx, args);
		derx0 = (that->*der)(x, args);
		derx1 = (that->*der)(x + dx, args);

		A[i] = fx0;
		B[i] = derx0;
		D[i] = (2 * (fx0 - fx1) + (derx0 + derx1) * dx) / dx;
		C[i] = (fx1 - fx0 + (-derx0 - D[i]) * dx);
	}
}

double Mesh::query(double x) {
	if(x <= xlow) return A[0];
	if(x >= xupp) x = xupp - FLT_EPSILON;
	int i = (int) ((x - xlow) / delta);
	double dx = x - xlow - delta * i;
	return (A[i] + dx * (B[i] + dx * (C[i] + dx * D[i]) * inv_sqr_delta));
}

double Mesh::query_derivative(double x) {
	if(x < xlow) return B[0];
	if(x >= xupp) x = xupp - FLT_EPSILON;
	int i = (int) ((x - xlow) / delta);
	double dx = x - xlow - delta * i;
	return (B[i] + (2 * dx * C[i] + 3 * dx * dx * D[i]) * inv_sqr_delta);
}

LookupTable::LookupTable(std::string filename, int points) : _lt_points(points) {
	std::ifstream lt_file(filename.c_str(), std::ios::in);
	if(!lt_file.good()) {
		std::cerr << "Can't read lookup file '" << filename << "'. Aborting" << std::endl;
		exit(1);
	}

	base_function data;
	std::string line;

	int i = 0;
	bool stop = false;
	while(!stop) {
		std::getline(lt_file, line);
		std::vector<std::string> spl = split(line, ' ');
		if(spl.size() != 3 || lt_file.eof()) stop = true;
		else {
			data.x[i] = atof(spl[0].c_str());
			data.fx[i] = atof(spl[1].c_str());
			data.dfx[i] = atof(spl[2].c_str());

			if(i > 0 && data.x[i] <= data.x[i - 1]) {
				std::cerr << "The x values of the lookup table should be monotonically increasing (found " << data.x[i] << " <= " << data.x[i - 1];
				exit(1);
			}
			i++;
		}
	}
	lt_file.close();

	double lowlimit = data.x[0];
	double uplimit = data.x[i - 1];

	_lookup_table.build(this, &LookupTable::_fx, &LookupTable::_dfx, (void *) (&data), _lt_points, lowlimit, uplimit);
}

LookupTable::~LookupTable() {

}

void LookupTable::potential(Bead &p, Bead &q, double shift_by) {
	double r = q.x - p.x + shift_by;

	double energy = _lookup_table.query(r);

	if(r < _lookup_table.xlow) {
		std::cout << "Exceeded the lower bound (" << r << " < " << _lookup_table.xlow << ")";
	}

	double force = -_lookup_table.query_derivative(r);
	p.force -= force;
	q.force += force;
}

double LookupTable::_linear_interpolation(double x, std::vector<double> &x_data, std::vector<double> fx_data, int points) {
	int ind = -1;
	for(int i = 0; i < points && ind == -1; i++)
		if(x_data[i] > x) ind = i;

	double val;
	if(ind == -1) {
		int last = points - 1;
		double slope = (fx_data[last] - fx_data[last - 1]) / (x_data[last] - x_data[last - 1]);
		val = fx_data[last] + slope * (x - x_data[last]);
	}
	else {
		double slope = (fx_data[ind] - fx_data[ind - 1]) / (x_data[ind] - x_data[ind - 1]);
		val = fx_data[ind - 1] + slope * (x - x_data[ind - 1]);
	}

	return val;
}

double LookupTable::_fx(double x, void *par) {
	base_function *data = (base_function *) par;
	return _linear_interpolation(x, data->x, data->fx, data->points);
}

double LookupTable::_dfx(double x, void *par) {
	base_function *data = (base_function *) par;
	return _linear_interpolation(x, data->x, data->dfx, data->points);
}
