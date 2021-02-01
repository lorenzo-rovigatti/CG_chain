#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>

#include "utils/strings.h"
#include "utils/Input.h"
#include "utils/thermostat.h"

#define SQR(x) ((x) * (x))
#define BEAD_SIZE 3

struct Chain;
struct Atom;
struct Input;

Chain build_from_topology_file(Input &input);
void LJ(Atom &p, Atom &q, double shift_by);

// a single atom
struct Atom {
	int bead_id;
	double m;
	double eps;
	double x;
	double v = 0.;
	double f_inter = 0.;
	double f_intra = 0.;
	double E = 0.;

	Atom(int bead, double nm, double neps, double nx) : bead_id(bead), m(nm), eps(neps), x(nx) {

	}

	double accel() {
		return (f_inter + f_intra) / m;
	}
};

// a bead (a collection of nearby BEAD_SIZE atoms)
struct Bead {
	std::array<Atom *, BEAD_SIZE> atoms;

	double R() {
		double mass = 0.;
		double pos = 0.;
		for(auto atom : atoms) {
			pos += atom->m * atom->x;
			mass += atom->m;
		}

		return pos / mass;
	}

	double F() {
		return atoms[BEAD_SIZE - 1]->f_inter;
	}
};

// the whole chain
struct Chain {
	std::vector<Atom> atoms;
	std::vector<Bead> beads;
	// a constant used to generate the initial configuration and to set the length of the chain
	static constexpr double bond_length = std::pow(2., 1. / 6.);

	int N() {
		return atoms.size();
	}

	double L() {
		return N() * bond_length;
	}
};

// this is a sampler which averages the simulation data and then prints out the final result
struct EffectiveForce {
	double min, max, bin_size;
	std::vector<double> data;
	std::vector<long long int> counter;
	
	EffectiveForce(double nmin, double nmax, int bins) :
		min(nmin),
		max(nmax),
		data(bins),
		counter(bins) {
		bin_size = (nmax - nmin) / (double) bins;
	}

	// accumulate the force
	void add(double R, double value) {
		if(R < min || R > max) {
			return;
		}

		int bin = (int) ((R - min) / bin_size);
		counter[bin]++;
		data[bin] += value;
	}

	// print the force as a function of the distance to stdout
	void print() {
		for(unsigned int i = 0; i < data.size(); i++) {
			if(counter[i] > 0) {
				double R = (i + 0.5) * bin_size + min;
				double F = data[i] / (double) counter[i];
				// we scale these values by the bond length to match Marco et al's units of measurements
				R /= Chain::bond_length;
				F *= Chain::bond_length;
				std::cout << R << " " << F << std::endl;
			}
		}
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
	EffectiveForce F_eff(2.5, 5., 100);
	int newtonian_steps = atoi(input.get("newtonian_steps"));
	double dt = atof(input.get("dt"));
	double diff_coeff = atof(input.get("diff_coeff"));
	long long int steps = atoll(input.get("steps"));
	long long int equilibration_steps = atoll(input.get("equilibration_steps"));
	long long int print_every = atoll(input.get("print_every"));
	long long int sample_every = atoll(input.get("sample_every"));

	double thermostat_pt = (2. * T * newtonian_steps * dt) / (T * newtonian_steps * dt  + 2 * diff_coeff);
	BrownianThermostat thermostat(T, thermostat_pt);

	// extract velocities from the correct maxwell boltzman distribution
	thermostat.random_velocities(chain.atoms);

	std::cerr << "INFO: thermostat pt: " << thermostat_pt << std::endl;

	std::ofstream energy_file("energy.dat");

	// run the simulation
	for(int i = 0; i < steps; i++) {
		// first integration step
		for(auto &atom : chain.atoms) {
			atom.v += atom.accel() * dt * 0.5;
			atom.x += atom.v * dt;
			atom.f_inter = atom.f_intra = atom.E = 0.;
		}

		// force calculation
		for(auto atom = chain.atoms.begin(); atom != chain.atoms.end(); atom++) {
			auto next = atom + 1;
			double shift_by = 0;
			if(next == chain.atoms.end()) {
				next = chain.atoms.begin();
				shift_by = chain.L();
			}
			LJ(*atom, *next, shift_by);
			
		}

		// second integration step
		for(auto &atom : chain.atoms) {
			atom.v += atom.accel() * dt * 0.5;
		}

		// we apply the thermostat once every newtonian_steps steps
		if(i % newtonian_steps == 0) {
			thermostat.apply(chain.atoms);
		}

		// print the potential, kinetic and total energy per particle
		if(i % print_every == 0) {
			// compute the centre of mass velocity
			double v_com = 0.;
			double m_com = 0.;
			for(auto &atom : chain.atoms) {
				v_com += atom.m * atom.v;
				m_com += atom.m;
				
			}
			v_com /= m_com;
			
			double pot_energy = 0.;
			double kin_energy = 0.;
			for(auto &atom : chain.atoms) {
				// remove the com velocity
				atom.v -= v_com;
				// we divide the potential energy of each atom by two because we assign the same energy to each particle in an interacting pair
				pot_energy += atom.E / 2.;
				kin_energy += 0.5 * atom.m * SQR(atom.v);
			}
			pot_energy /= chain.N();
			kin_energy /= chain.N();
			double energy_per_particle = pot_energy + kin_energy;
			energy_file << i * dt << " " << pot_energy << " " << kin_energy << " " << energy_per_particle << std::endl;
		}

		// sample the inter-bead force
		if(i > equilibration_steps && i % sample_every == 0) {
			for(auto bead = chain.beads.begin(); bead != chain.beads.end(); bead++) {
				auto next = bead + 1;
				double shift_by = 0.;
				if(next == chain.beads.end()) {
					next = chain.beads.begin();
					shift_by = chain.L();
				}
				double R = next->R() - bead->R() + shift_by;

				// the minus sign comes from the fact that we take the force acting on the bead that is on the left: if the force
				// acting on it is negative then the bead will be pushed towards larger R (i.e. itis repelled by the other bead and
				// hence its "force" should be positive)
				F_eff.add(R, -bead->F());
			}
		}
	}

	energy_file.close();

	// print the force as a function of the inter-bead distance
	F_eff.print();
	
	return 0;
}

void LJ(Atom &p, Atom &q, double shift_by) {
	double r = q.x - p.x + shift_by;
	double eps = std::sqrt(p.eps * q.eps);
	double ir6 = 1. / (SQR(r) * SQR(r) * SQR(r));
	double energy = 4. * eps * (SQR(ir6) - ir6);
	double force = (24. * eps * (ir6 - 2. * SQR(ir6))) / r;

	p.E += energy;
	q.E += energy;

	if(p.bead_id == q.bead_id) {
		p.f_intra += force;
		q.f_intra -= force;
	}
	else {
		p.f_inter += force;
		q.f_inter -= force;
	}
}

Chain build_from_topology_file(Input &input) {
	const char *filename = input.get("topology_file");
	std::ifstream inp(filename);

	if(!inp.good()) {
		std::cerr << "Invalid topology file " << filename << std::endl;
		exit(1);
	}

	Chain chain;

	int N_atoms;
	inp >> N_atoms;

	if((N_atoms % BEAD_SIZE) != 0) {
		std::cerr << "The number of atoms should be a multiple of the bead size " << BEAD_SIZE << std::endl;
		exit(1);
	}

	chain.atoms.reserve(N_atoms);
	chain.beads.resize(N_atoms / BEAD_SIZE);
	int curr_bead_id = 0;
	for(int i = 0; i < N_atoms; i++) {
		if(i > 0 && i % (BEAD_SIZE) == 0) {
			curr_bead_id++;
		}
		double x = chain.bond_length * i;
		std::string species;
		inp >> species;
		double mass, eps;
		std::tie<double, double>(mass, eps) = input.get_mass_and_epsilon(species);

		Atom atom(curr_bead_id, mass, eps, x);
		atom.v = atom.f_inter = atom.f_intra = 0.;
		chain.atoms.push_back(atom);
		chain.beads[curr_bead_id].atoms[i % BEAD_SIZE] = &(chain.atoms.data()[i]);
	}

	inp.close();

	return chain;
}

