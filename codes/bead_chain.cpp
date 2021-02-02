#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>

#include "CG/definitions.h"
#include "CG/LookupTable.h"
#include "utils/strings.h"
#include "utils/Input.h"
#include "utils/thermostat.h"

struct Chain;
struct Input;

Chain build_from_topology_file(Input &input);

// the whole chain
struct Chain {
	int _bead_size;
	std::vector<CGBead> beads;
	// a constant used to generate the initial configuration and to set the length of the chain
	static constexpr double bond_length = std::pow(2., 1. / 6.);

	Chain(int bead_size) : _bead_size(bead_size) {

	}

	int N() {
		return beads.size();
	}

	double bead_bond_length() {
		return _bead_size * bond_length;
	}

	double L() {
		return N() * bead_bond_length();
	}
};

int main(int argc, char *argv[]) {
	if(argc < 2) {
		std::cerr << "Usage is " << argv[0] << " input_file" << std::endl;
		exit(1);
	}

	srand48(591823);

	// parse the input file
	Input input(argv[1]);

	// get the values we need from the input file
	double T = 1. / atof(input.get("beta"));
	Chain chain = build_from_topology_file(input);

	int newtonian_steps = atoi(input.get("newtonian_steps"));
	double dt = atof(input.get("dt"));
	double diff_coeff = atof(input.get("diff_coeff"));
	long long int steps = atoll(input.get("steps"));
	long long int equilibration_steps = atoll(input.get("equilibration_steps"));
	long long int print_every = atoll(input.get("print_every"));
	long long int sample_every = atoll(input.get("sample_every"));
	long long int print_bead_conf_every = atoll(input.get("print_bead_conf_every", "0"));
	std::string trajectory_file(input.get("bead_trajectory_file", "trajectory.dat"));
	if(print_bead_conf_every) {
		std::ofstream trajectory(trajectory_file.c_str());
		trajectory.close();
	}

	double thermostat_pt = (2. * T * newtonian_steps * dt) / (T * newtonian_steps * dt  + 2 * diff_coeff);
	BrownianThermostat thermostat(T, thermostat_pt);

	LookupTable lt(input.get("lookup_table_file"), 100);

	// extract velocities from the correct maxwell boltzman distribution
	thermostat.random_velocities(chain.beads);

	std::cerr << "INFO: thermostat pt: " << thermostat_pt << std::endl;

	std::ofstream energy_file("energy.dat");

	// run the simulation
	for(int i = 0; i < steps; i++) {
		// first integration step
		for(auto &bead : chain.beads) {
			bead.v += bead.accel() * dt * 0.5;
			bead.x += bead.v * dt;
			bead.force = bead.E = 0.;
		}

		// force calculation
		for(auto bead = chain.beads.begin(); bead != chain.beads.end(); bead++) {
			auto next = bead + 1;
			double shift_by = 0;
			if(next == chain.beads.end()) {
				next = chain.beads.begin();
				shift_by = chain.L();
			}
			lt.potential(*bead, *next, shift_by);
			
		}

		// second integration step
		for(auto &bead : chain.beads) {
			bead.v += bead.accel() * dt * 0.5;
		}

		// we apply the thermostat once every newtonian_steps steps
		if(i % newtonian_steps == 0) {
			thermostat.apply(chain.beads);
		}

		// print the potential, kinetic and total energy per particle
		if(i % print_every == 0) {
			// compute the centre of mass velocity
			double v_com = 0.;
			double m_com = 0.;
			for(auto &bead : chain.beads) {
				v_com += bead.m * bead.v;
				m_com += bead.m;
				
			}
			v_com /= m_com;
			
			double pot_energy = 0.;
			double kin_energy = 0.;
			for(auto &bead : chain.beads) {
				// remove the com velocity
				bead.v -= v_com;
				// we divide the potential energy of each atom by two because we assign the same energy to each particle in an interacting pair
				pot_energy += bead.E / 2.;
				kin_energy += 0.5 * bead.m * SQR(bead.v);
			}
			pot_energy /= chain.N();
			kin_energy /= chain.N();
			double energy_per_particle = pot_energy + kin_energy;
			energy_file << i * dt << " " << pot_energy << " " << kin_energy << " " << energy_per_particle << std::endl;
		}


		if(i > equilibration_steps && print_bead_conf_every && (i % print_bead_conf_every == 0)) {
			std::ofstream trajectory(trajectory_file.c_str(), std::ios::app);

			trajectory << chain.beads.size() << " " << chain.L() << " " << i << std::endl;
			for(auto &bead : chain.beads) {
				trajectory << bead.x << std::endl;
			}

			trajectory.close();
		}
	}

	energy_file.close();

	return 0;
}

Chain build_from_topology_file(Input &input) {
	const char *filename = input.get("topology_file");
	std::ifstream inp(filename);

	if(!inp.good()) {
		std::cerr << "Invalid topology file " << filename << std::endl;
		exit(1);
	}

	double mass = atof(input.get("bead_mass"));
	int bead_size = atoi(input.get("bead_size"));

	Chain chain(bead_size);

	int N_beads;
	inp >> N_beads;

	chain.beads.reserve(N_beads);
	int curr_bead_id = 0;
	for(int i = 0; i < N_beads; i++) {
		double x = chain.bead_bond_length() * i;

		CGBead bead(mass, x);

		chain.beads.push_back(bead);
	}

	inp.close();

	return chain;
}

