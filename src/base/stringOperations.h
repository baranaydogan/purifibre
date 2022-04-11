#ifndef SRC_STRINGOPERATIONS_H_
#define SRC_STRINGOPERATIONS_H_

#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>

std::vector<std::string> splitString (const std::string &s, char delim);

std::string getFileExtension(std::string filePath);
std::string removeFileExtension(std::string filePath);

bool existsFile(const std::string& name);

void im(std::string m);

void msg_error(std::string m);

#endif
