#ifndef SRC_IMAGE_SF_IMAGE_H_
#define SRC_IMAGE_SF_IMAGE_H_

#include "../math/core.h"
#include "../math/sphericalFunctions.h"
#include "image.h"
#include "vector"

class SF_Image : public Image<float> {
    
public:

    SF_Image();
    SF_Image(std::string _filePath, bool _isEven);
    SF_Image(const char* _filePath, bool _isEven) {SF_Image(std::string(_filePath),_isEven);}
    ~SF_Image();

    bool    read();
    float   getSFval(float *p, float* tan) {return (*this)(p,SF::coordinate2index(tan));}
    void    smooth(float angle);

private:

    bool    isEven;
    int     sphericalDomainResolution;

};

#endif
