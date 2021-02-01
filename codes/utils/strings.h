/*
 * strings.h
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#ifndef CODES_UTILS_STRINGS_H_
#define CODES_UTILS_STRINGS_H_

#include <vector>
#include <string>

std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");

std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");

std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ");

std::vector<std::string> split(const std::string &s, char delim);

#endif /* CODES_UTILS_STRINGS_H_ */
