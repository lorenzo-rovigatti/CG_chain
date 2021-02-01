/*
 * definitions.h
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#ifndef CODES_CG_DEFINITIONS_H_
#define CODES_CG_DEFINITIONS_H_

// a single atom
struct Bead {
	double m;
	double x;
	double v = 0.;
	double force = 0.;
	double E = 0.;

	Bead(double nm, double nx) : m(nm), x(nx) {

	}

	double accel() {
		return force / m;
	}
};

#endif /* CODES_CG_DEFINITIONS_H_ */
