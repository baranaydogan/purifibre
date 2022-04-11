#include "stringOperations.h"

std::vector<std::string> splitString (const std::string &s, char delim) {
    
    std::vector<std::string> out;
    
    std::stringstream ss(s);
    
    std::string item;

    while (getline (ss, item, delim)) {
        out.push_back (item);
    }

    return out;
}


std::string getFileExtension(std::string filePath) {

    size_t lastDot = filePath.find_last_of(".");
    
    if (lastDot==std::string::npos)
        return "";
    
    std::string ext = filePath.substr(lastDot+1);
    
    std::string rem = filePath.substr(0,lastDot);
    size_t dotCheck = rem.find_last_of(".");
        
    if (dotCheck==std::string::npos)
        return ext;
    
    if (ext=="gz") {
        ext = rem.substr(dotCheck+1);
        return (ext + ".gz");
    } else
        return ext;
    
}

std::string removeFileExtension(std::string filePath) {

    size_t lastDot = filePath.find_last_of(".");
    
    if (lastDot==std::string::npos)
        return "";
    
    std::string ext = filePath.substr(lastDot+1);
    
    std::string rem = filePath.substr(0,lastDot);
    size_t dotCheck = rem.find_last_of(".");
        
    if (dotCheck==std::string::npos)
        return rem;
    
    if (ext=="gz") {
        return rem.substr(0,dotCheck);
    } else
        return rem;
    
}

// Returns true if file exists
bool existsFile(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}

void im(std::string m) {std::cout << m << std::endl << std::flush;}

void msg_error(std::string m) {std::cout << m << std::endl << std::flush;}
