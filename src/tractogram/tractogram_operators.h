#ifndef SRC_TRACTOGRAM_OPERATORS_H_
#define SRC_TRACTOGRAM_OPERATORS_H_

#include "tractogramReader.h"

// Compute bounding box of a tractogram
std::vector<float> getTractogramBBox(TractogramReader* tractogram);


#endif