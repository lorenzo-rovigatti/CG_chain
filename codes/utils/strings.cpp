/*
 * strings.cpp
 *
 *  Created on: Feb 1, 2021
 *      Author: lorenzo
 */

#include "strings.h"

#include <sstream>

std::string& ltrim(std::string& str, const std::string& chars) {
    str.erase(0, str.find_first_not_of(chars));
    return str;
}

std::string& rtrim(std::string& str, const std::string& chars) {
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}

std::string& trim(std::string& str, const std::string& chars) {
    return ltrim(rtrim(str, chars), chars);
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::string s_copy(s);
	trim(s_copy);
	std::vector<std::string> elems;
	std::stringstream ss(s_copy);
	std::string item;

	while(getline(ss, item, delim)) {
		trim(item);
		if(item.length() > 0) {
			elems.push_back(item);
		}
	}

	return elems;
}
