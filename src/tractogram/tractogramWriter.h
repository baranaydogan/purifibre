#ifndef SRC_TRACTOGRAMWRITER_H_
#define SRC_TRACTOGRAMWRITER_H_

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <float.h>
#include "../base/config.h"
#include "../base/dataTypeHandler.h"
#include "../base/vectorOperations.h"
#include "../math/core.h"
#include "tractogramReader.h"


typedef enum {
    STREAMLINE_OWNER,
    POINT_OWNER,
    OWNER_NOTSET
} TractogramOwnerType;

struct TractogramField {
    TractogramOwnerType     owner;
    std::string             name;
    DATATYPE                datatype;
    int                     dimension;
    std::vector<void*>      data;
};

bool writeTractogram(std::string fname,std::vector<std::vector<std::vector<float>>>& tractogram);
bool writeTractogram(std::string fname,std::vector<std::vector<std::vector<float>>>& tractogram,std::vector<TractogramField>& fields);

bool writeTractogram(std::string out_fname,TractogramReader* tractogram,std::vector<size_t>& idx);
bool writeTractogram(std::string out_fname,TractogramReader* tractogram);
bool writeTractogram(std::string out_fname,std::string inp_fname,std::vector<size_t>& idx);
bool writeTractogram(std::string out_kept_fname,std::string out_rmvd_fname,std::string inp_fname,std::vector<size_t>& idx);

#endif
