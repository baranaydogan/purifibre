#ifndef SRC_MATH_SPHERICALFUNCTIONS_H_
#define SRC_MATH_SPHERICALFUNCTIONS_H_

#include "core.h"

namespace SF {

extern std::vector<std::vector<float>>                      sfCoords;
extern std::vector<std::vector<std::tuple<int,float> >>     sfNeighbors; // Contains indices of neighbors for each coordinate in the order of increasing distance 
extern bool                                                 sfIsEven;

void    init(bool _sfIsEven, int _sfDim);
void    clean();
int64_t coordinate2index(float* coord);
void    coordinateNeighbors(std::vector<std::tuple<int,float> >& neighbors, float* coord, float distThresh);

}

#endif
