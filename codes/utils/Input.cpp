/*
 * Input.cpp
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#include "Input.h"

#include "strings.h"

#include <fstream>

Input::Input(char *filename) {
	std::ifstream inp(filename);

	if(!inp.good()) {
		std::cerr << "Invalid input file " << filename << std::endl;
		exit(1);
	}

	while(inp.good()) {
		std::string line;
		std::getline(inp, line);
		line = trim(line);
		// skip empty lines and lines starting with #
		if(line.size() == 0 || line[0] == '#') {
			continue;
		}
		auto spl = split(line, '=');
		if(spl.size() != 2) {
			std::cerr << "Invalid line '" << line << "'" << std::endl;
			exit(1);
		}

		input_map[spl[0]] = spl[1];
	}
}

// get a key from the stored map, exiting if the key is not found
const char *Input::get(std::string key, std::string default_value) {
	if(input_map.find(key) == input_map.end()) {
		if(default_value.length() > 0) {
			return default_value.c_str();
		}
		else {
			std::cerr << "Key '" << key << "' not found in the input file" << std::endl;
			exit(1);
		}
	}

	return input_map[key].c_str();
}

// keys that specify a species' mass and epsilon are treated differently. This method returns a pair
// containing the mass and epsilon values associated to the specific key (which should be the name of the species)
std::pair<double, double> Input::get_mass_and_epsilon(std::string key) {
	std::string value = get(key);
	auto spl = split(value, ',');
	double mass = atof(spl[0].c_str());
	double epsilon = atof(spl[1].c_str());

	return std::pair<double, double>({mass, epsilon});
}

