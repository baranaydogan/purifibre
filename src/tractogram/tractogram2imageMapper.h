#ifndef SRC_TRACTOGRAM2IIMAGEMAPPER_H_
#define SRC_TRACTOGRAM2IIMAGEMAPPER_H_

#include "../math/core.h"
#include "tractogramReader.h"
#include "../image/image.h"

template<typename T>
class Tractogram2ImageMapper {

public:

    Tractogram2ImageMapper(TractogramReader* _tractogram, Image<T>* _img);
    Tractogram2ImageMapper(TractogramReader* _tractogram, Image<T>* _img, bool allocateGrid);
    
    ~Tractogram2ImageMapper();
    
    void run (
        std::function<void(Tractogram2ImageMapper<T>* tim, uint16_t _threadNo, int* _gridPos, Segment& _seg)> processor_f,
        std::function<void(Tractogram2ImageMapper<T>* tim)> outputCompiler_f
    );

    bool processStreamline(int _streamlineId, uint16_t _threadNo, std::function<void(Tractogram2ImageMapper<T>* tim, uint16_t _threadNo, int* _gridPos, Segment& _seg)> f );
    
    void setMapOnce(bool val) {mapOnce=val;}
    bool setMask(Image<int>* maskImg);
    bool setMask(Image<int>* maskImg, int selectedLabel);
    void setMask(bool*** _mask) {mask = _mask;}
    void anisotropicSmoothing(std::tuple<float,int> _smoothing) {smoothing = _smoothing;}
    
    std::vector<std::vector<std::vector<std::vector<void*>>>> grid; // Grid holds a void* for each [threadNo][i][j][k]
    Image<T>*         img;

private:
    
    TractogramReader* tractogram;

    bool***                 mask;
    bool                    mapOnce;
    bool                    gridAllocated;
    std::tuple<float,int>   smoothing;
    
};

#endif