# CG_chain

## Compilation

Download and compile the codes with

	$ git clone https://github.com/lorenzo-rovigatti/CG_chain.git
	$ cd CG_chain
	$ mkdir build
	$ cd build
	$ cmake ../codes
	$ make
    
Two executables will be created: `Veff_Chain` (for atomistic simulations) and `bead_chain` (for CG simulations). Simulations can be run by passing an input file to an executable:

	$ ./Veff_chain input > force.dat # start an atomistic simulation and print the resulting effective force to the force.dat file 
	$ ./bead_chain input # start a CG simulation
    
Both executables will produce an `energy.dat` file containing three four columns: the time, potential, kinetic and total energy. The executables can optionally print also the configuration files (see below for details). 

At the end of the simulation `Veff_chain` will print the effective force. Redirect its output to file to save it.
    
## Input files

Options common for both types of simulations are:

* `beta`: inverse temperature. Set it to `1`.
* `topology_file`: a file containing the topology of the chains we want to simulate. See below.
* `steps`: number of steps that the simulation should run for.
* `equilibration_steps`: number of steps during which no sampling (configuration printing / effective-force computing) will occur.
* `print_every`: number of steps every which the energy should be printed (reasonable values are 1000 or 10000).
* `dt`: integration time steps. `0.001` is a good value.
* `diff_coeff`: set it to `0.1` or `0.5`. Used by the thermostat and not so important for structural values.
* `newtonian_steps`: set it to `53` or something like that, another thermostat parameter.
* `print_bead_conf_every`: frequency with which configurations should be printed. Defaults to 0 (i.e. do not print configurations).
* `bead_conf_prefix`: configurations will be printed to files having this prefix, followed by the current time step. It can (and should, if you want my advise) contain folders (`confs/bead_conf_`, for instance).

### Atomistic simulations

Specific options:

* `sample_every`: number of steps every which the inter-bead force will be sampled. Use `100` or something like that.
* `force_rescale`: if not 0 the effective force and will be rescaled by the bond length. Useful to compare with the Julia code. It defaults to 0 (i.e. do not rescale).

In the input files you also need to specify the mass and interaction constant of the different species. This is done by setting special keys that contain two comma-separated values. If you want to simulate two atoms of mass 1 and 10 with interaction constants 10 and 0.1 you would add these two options to the input file:

	A = 1, 10
	B = 10, 0.1

Note that the name of the options (here `A` and `B`) need to be the same used in the topology file (see below)

### CG simulations

Specific options:

* `bead_mass`: total mass of a bead. It should not be important if you are interested in structural quantities.
* `bead_size`: number of atoms in a bead. We usually coarse grain 3 atoms, should you should set it to 3.
* `lookup_table_file`: a file containing three columns: `x U(x) dU/dx`, where `x` is the distance between two beads, `U(x)` the corresponding energy and `dU/dx` its derivative. 

## Topology files

### Atomistic simulations

The first line should contain the number of atoms `N` you want to simulate, and the `N` lines that follow should each contain the type of atom they refer to. A chain made of 6 atoms in a `ABAABA' sequence can be simulated with the following topology file:

	6
	A
	B
	A
	A
	B
	A

### CG simulations

The topology should contain only one value: the number of beads in the chain.
