/*
 * thermostat.h
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#ifndef CODES_UTILS_THERMOSTAT_H_
#define CODES_UTILS_THERMOSTAT_H_

// get a random number from a gaussian with zero mean and unity variance
double gaussian() {
	static unsigned int isNextG = 0;
	static double nextG;
	double toRet;
	double u, v, w;

	if(isNextG) {
		isNextG = 0;
		return nextG;
	}

	w = 2.;
	while(w >= 1.0) {
		u = 2. * drand48() - 1.0;
		v = 2. * drand48() - 1.0;
		w = u * u + v * v;
	}

	w = std::sqrt((-2. * std::log(w)) / w);
	toRet = u * w;
	nextG = v * w;
	isNextG = 1;

	return toRet;
}

class BrownianThermostat {
public:
	BrownianThermostat(double T, double pt) {
		_T = T;
		_pt = pt;
	}

	template <typename particle_vector>
	void random_velocities(particle_vector &v) {
		double finite_size_factor = v.size() / (v.size() - 1.);

		for(auto &p : v) {
			double factor = std::sqrt(_T / p.m * finite_size_factor);
			p.v = gaussian() * factor;
		}
	}

	template <typename particle_vector>
	void apply(particle_vector &v) {
		double finite_size_factor = v.size() / (v.size() - 1.);

		for(auto &p : v) {
			// there is a certain chance (pt) that an atom undergoes a brownian collision so that its velocity is extracted anew from a maxwellian
			if(drand48() < _pt) {
				double factor = std::sqrt(_T / p.m * finite_size_factor);
				p.v = gaussian() * factor;
			}
		}
	}

private:
	double _T;
	double _pt;
};

#endif /* CODES_UTILS_THERMOSTAT_H_ */
