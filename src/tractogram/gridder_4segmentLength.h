#ifndef SRC_GRIDDER_4LENGTH_H_
#define SRC_GRIDDER_4LENGTH_H_

#include <stdio.h>
#include <map>
#include "../math/sphericalFunctions.h"
#include "tractogram2imageMapper.h"

using namespace SF;

// Regular output
template<typename T>
inline void processor_4segmentLength(Tractogram2ImageMapper<T>* tim, uint16_t threadNo, int* gridPos, Segment& seg) {
    if (tim->grid[threadNo][gridPos[0]][gridPos[1]][gridPos[2]]==NULL) {                    // Allocate memory here
        float* segmentLength                                    = new float;
        *segmentLength                                          = 0;
        tim->grid[threadNo][gridPos[0]][gridPos[1]][gridPos[2]] = ((void*)segmentLength);   // Each thread handles a single float for each voxel in the output
    }
    *((float*)(tim->grid[threadNo][gridPos[0]][gridPos[1]][gridPos[2]])) += seg.length;     // Add to the segment length
}

template<typename T>
inline void outputCompiler_4segmentLength(Tractogram2ImageMapper<T>* tim) {

    tim->img->allocData();

    auto genOut = [&](MTTASK task)->void {
        int64_t ind = 0;
        for (int i = 0; i < tim->img->imgDims[0]; i++)
            for (int j = 0; j < tim->img->imgDims[1]; j++)
                for (int t = 0; t < MT::maxNumberOfThreads; t++)
                    if (tim->grid[t][i][j][task.no]!=NULL) {
                        ind = tim->img->sub2ind(i,j,task.no);
                        tim->img->data[ind] += *((float*)(tim->grid[t][i][j][task.no]));
                        delete (float*)(tim->grid[t][i][j][task.no]);                       // Deallocate memory here
                    }
    };
    MT::MTRUN(tim->img->imgDims[2], genOut); // Because addition is not thread-safe, we cannot write on the same voxel with different threads. So we multi-thread along z-dimension.

}

// Spherical output
template<typename T>
inline void processor_4segmentLength_sf(Tractogram2ImageMapper<T>* tim, uint16_t threadNo, int* gridPos, Segment& seg) {
    if (tim->grid[threadNo][gridPos[0]][gridPos[1]][gridPos[2]]==NULL) {
        std::map<int64_t,float>* dirLength = new std::map<int64_t,float>();
        tim->grid[threadNo][gridPos[0]][gridPos[1]][gridPos[2]] = ((void*)dirLength);   // Allocate memory here
    }
    (*((std::map<int64_t,float>*)(tim->grid[threadNo][gridPos[0]][gridPos[1]][gridPos[2]])))[SF::coordinate2index(seg.dir)] += seg.length;
}

template<typename T>
inline void outputCompiler_4segmentLength_sf(Tractogram2ImageMapper<T>* tim) {

    tim->img->allocData();

    auto genOut = [&](MTTASK task)->void {
        int64_t ind = 0;
        for (int i = 0; i < tim->img->imgDims[0]; i++)
            for (int j = 0; j < tim->img->imgDims[1]; j++)
                for (int t = 0; t < MT::maxNumberOfThreads; t++)
                    if (tim->grid[t][i][j][task.no]!=NULL) {
                        
                        for (auto lenSf : (*((std::map<int64_t,float>*)(tim->grid[t][i][j][task.no])))) {
                            ind = tim->img->sub2ind(i,j,task.no,lenSf.first);
                            tim->img->data[ind] += lenSf.second;
                        }

                        delete (std::map<int64_t,float>*)(tim->grid[t][i][j][task.no]); // Deallocate memory here
                    }
    };

     // Because addition is not thread-safe, we cannot write on the same voxel with different threads. So we multi-thread along z-dimension.
    if (!QUITE)
        MT::MTRUN(tim->img->imgDims[2], "Compiling spherical function", genOut);
    else
        MT::MTRUN(tim->img->imgDims[2], genOut);

}

#endif
