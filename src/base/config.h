#ifndef SRC_CONFIG_H_
#define SRC_CONFIG_H_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>


#include "CLI11.hpp"
#include "byteSwapper.h"
#include "multithreader.h"
#include "stringOperations.h"
#include "dataTypeHandler.h"

using namespace MT;

extern std::string SGNTR;
extern std::string FOOTER;
extern bool        QUITE;

extern void INIT();

typedef enum {
    UNKNOWNSPACEUNIT,
    METER,
    MM,
    MICRON
} SPACEUNIT;

typedef enum {
    UNKNOWNTIMEUNIT,
    SEC,
    MSEC,
    USEC,
    HZ,
    PPM,
    RADS
} TIMEUNIT;


#endif
