/*
 * Input.h
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#ifndef CODES_UTILS_INPUT_H_
#define CODES_UTILS_INPUT_H_

#include <map>
#include <iostream>

// a very simple input parser
struct Input {
	std::map<std::string, std::string> input_map;

	// parse the input file
	Input(char *filename);

	// get a key from the stored map, exiting if the key is not found
	const char *get(std::string key, std::string default_value="");

	// keys that specify a species' mass and epsilon are treated differently. This method returns a pair
	// containing the mass and epsilon values associated to the specific key (which should be the name of the species)
	std::pair<double, double> get_mass_and_epsilon(std::string key);
};

#endif /* CODES_UTILS_INPUT_H_ */
