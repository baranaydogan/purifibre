#ifndef SRC_IMAGE_CPP_
#define SRC_IMAGE_CPP_

#include "image.h"

// Explicit instantiations
template class Image<bool>;
template class Image<uint8_t>;
template class Image<int8_t>;
template class Image<uint16_t>;
template class Image<int16_t>;
template class Image<uint32_t>;
template class Image<int32_t>;
template class Image<uint64_t>;
template class Image<int64_t>;
template class Image<float>;
template class Image<double>;
template class Image<long double>;

// TODO: Implement converters for complex data types in image_reader.cpp
// template class Image<std::complex<float>>;
// template class Image<std::complex<double>>;
// template class Image<std::complex<long double>>;



template<typename T>
void Image<T>::init()
{
    filePath 		= "";
    fileExtension   = "";
    description     = "";
    numberOfDimensions = 0;
    memset(imgDims,0,7*sizeof(int64_t));
    memset(pixDims,0,7*sizeof(float));
    dataType        = TYPEIDS[typeid(T)];
    inputDataType   = dataType;
    dataScaler      = 1;
    dataOffset      = 0;
    spaceUnit       = UNKNOWNSPACEUNIT;
    timeUnit        = UNKNOWNTIMEUNIT;
    headerIsRead    = false;
    data            = NULL;
    outsideVal      = 0;
    interpMethod    = LINEAR;
    setInterpolationMethod(LINEAR);
    for (auto i=0; i<7; i++) indexOrder[i]=i;
}


template<typename T>
Image<T>::Image() {
    init();
}

template<typename T>
Image<T>::Image(std::string _filePath) {
    init();
    setFilePath(_filePath);
    headerIsRead    = readHeader();
}

template<typename T>
Image<T>::Image(const char* _filePath) {
    init();
    setFilePath(std::string(_filePath));
    headerIsRead    = readHeader();
}

template<typename T>
Image<T>::Image(std::string _filePath, int* _indexOrder) {
    init();
    setFilePath(_filePath);
    for (auto i=0; i<7; i++) indexOrder[i]=_indexOrder[i];
    headerIsRead    = readHeader();
}

template<typename T>
Image<T>::Image(const char* _filePath, int* _indexOrder) {
    init();
    setFilePath(std::string(_filePath));
    for (auto i=0; i<7; i++) indexOrder[i]=_indexOrder[i];
    headerIsRead    = readHeader();
}

template<typename T>
Image<T>::~Image() {
    if (data!=NULL) {
        delete[] data;
        data = NULL;
    };
}


template<typename T>
Image<T>::Image(const Image &img) {
    
    filePath            = img.filePath;
    fileExtension       = img.fileExtension;
    description         = img.description;

    spaceUnit           = img.spaceUnit;
    timeUnit            = img.timeUnit;
    voxCnt              = img.voxCnt;
    valCnt              = img.valCnt;
    numberOfDimensions  = img.numberOfDimensions;
    memcpy(imgDims,img.imgDims,7*sizeof(int64_t));
    memcpy(pixDims,img.pixDims,7*sizeof(float));
    smallestPixDim      = img.smallestPixDim;
    memcpy(xyz2ijk, img.xyz2ijk, 12*sizeof(float));
    memcpy(ijk2xyz, img.ijk2xyz, 12*sizeof(float));
    dataType            = img.dataType;
    memcpy(rastkr2ras, img.rastkr2ras, 16*sizeof(float));
    interpMethod        = img.interpMethod;
    outsideVal          = img.outsideVal;
    headerIsRead        = img.headerIsRead;
    memcpy(s2i, img.s2i, 7*sizeof(int64_t));
    dataScaler          = img.dataScaler;
    dataOffset          = img.dataOffset;
    memcpy(indexOrder,img.indexOrder,7*sizeof(int));
    inputDataType       = img.inputDataType;
    setInterpolationMethod(interpMethod);
    data                = NULL;
}

template<typename T>
void Image<T>::create(int ndims, int64_t* _imgDims, float* _pixDims, float _ijk2xyz[][4], bool allocateData)
{

    if (ndims>7)
        msg_error("Number of dimensions cannot be bigger than 7"); 

    numberOfDimensions = ndims;

    for (int i=0; i<7; i++) {
        imgDims[i] = 1;
        pixDims[i] = 1;
    }

    for (int i=0; i<ndims; i++) {
        imgDims[i] = _imgDims[i];
        pixDims[i] = _pixDims[i];
    }

    if (_ijk2xyz != NULL) {

        // Compute xyz2ijk
        mat44 ijk2xyz_m44;
        for (int i=0; i<3; i++)
            for (int j=0; j<4; j++) {
                ijk2xyz[i][j]       = _ijk2xyz[i][j];
                ijk2xyz_m44.m[i][j] =  ijk2xyz[i][j];
            }
        ijk2xyz_m44.m[3][0] = ijk2xyz_m44.m[3][1] = ijk2xyz_m44.m[3][2] = 0;
        ijk2xyz_m44.m[3][3] = 1;

        mat44 xyz2ijk_m44;
        xyz2ijk_m44 = nifti_mat44_inverse(ijk2xyz_m44);
        for (int i=0; i<3; i++)
            for (int j=0; j<4; j++) {
                xyz2ijk[i][j] = xyz2ijk_m44.m[i][j];
            }
    } else {
        for (int i=0; i<3; i++) {
            for (int j=0; j<4; j++) {
                xyz2ijk[i][j] = 0;
                ijk2xyz[i][j] = 0;
            }
            xyz2ijk[i][i] = 1;
            ijk2xyz[i][i] = 1;
        }
    }

    parseHeader();

    // Allocate data array
    if (allocateData)
        data = new T[voxCnt*valCnt]();
    else
        data = NULL;

}

template<typename T>
void Image<T>::createFromTemplate(Image* img, bool allocateData)
{
    create(img->numberOfDimensions, img->imgDims, img->pixDims, img->ijk2xyz, allocateData);
}

template<typename T>
void Image<T>::createFromTemplate(Image* img, float pixDimScaler, bool allocateData)
{

    if (pixDimScaler==1.0) {
        create(img->numberOfDimensions, img->imgDims, img->pixDims, img->ijk2xyz, allocateData);
        return;
    }

    float range[3];
    float scale[3];

    memcpy(imgDims, img->imgDims, 7*sizeof(int64_t));
    memcpy(pixDims, img->pixDims, 7*sizeof(float));

    for (int i=0; i<3; i++) {
        range[i]   = img->pixDims[i] * img->imgDims[i];
        imgDims[i] = range[i] / (img->pixDims[i] / pixDimScaler);
        pixDims[i] = range[i] / imgDims[i];
        scale[i]   = pixDims[i] / img->pixDims[i];
    }

    memcpy(ijk2xyz, img->ijk2xyz, 12*sizeof(float));
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            ijk2xyz[i][j] *= scale[i];

    // Compute xyz2ijk
    mat44 ijk2xyz_m44;
    for (int i=0; i<3; i++)
        for (int j=0; j<4; j++)
            ijk2xyz_m44.m[i][j] = ijk2xyz[i][j];

    mat44 xyz2ijk_m44;
    xyz2ijk_m44 = nifti_mat44_inverse(ijk2xyz_m44);
    for (int i=0; i<3; i++)
        for (int j=0; j<4; j++)
            xyz2ijk[i][j] = xyz2ijk_m44.m[i][j];

    numberOfDimensions = img->numberOfDimensions;

    parseHeader();

    // Allocate data array
    if (allocateData)
        data = new T[voxCnt*valCnt]();
    else
        data = NULL;
}

template<typename T>
void Image<T>::createFromBoundingBox(int ndim, std::vector<float> bb, bool allocateData) {

    createFromBoundingBox(ndim, bb, 1, allocateData);
    
}

template<typename T>
void Image<T>::createFromBoundingBox(int ndim, std::vector<float> bb, float _pixDim, bool allocateData) {

    std::vector<float> pd;
    pd.push_back(_pixDim);

    createFromBoundingBox(ndim, bb, pd, allocateData);
    
}

template<typename T>
void Image<T>::createFromBoundingBox(int ndim, std::vector<float> bb, std::vector<float> _pixDims, bool allocateData) {

    int64_t  imgD[7] = {1,1,1,1,1,1,1};
    float    pixD[7] = {1,1,1,1,1,1,1};
    float    _ijk2xyz[3][4];
    float    shift[3];

    if (_pixDims.size()==1) {
       for (int i=0; i<3; i++)
            pixD[i] = _pixDims[0];
    } else if (!_pixDims.empty()) {
        for (int i=0; i<ndim; i++)
            pixD[i] = _pixDims[i];
    }


    for (int i=0; i<ndim; i++) {

        float length  = bb[i*2+1]-bb[i*2];
        float center  = length*0.5f;

        imgD[i]       = int64_t(length/pixD[i]);
        
        if ((imgD[i]*pixD[i]) < length) 
            imgD[i]++;

        float digiLen = imgD[i]*pixD[i];
        float digiCen = digiLen*0.5f;

        if (i<3) 
            shift[i] = bb[i*2] + (center - digiCen) + pixD[i]*0.5f;
        
        if (imgD[i]<0) 
            imgD[i] = 1;
    }

    

    for (int i=0; i<3; i++) {
        for (int j=0; j<4; j++) {
            _ijk2xyz[i][j] = 0;
        }
        _ijk2xyz[i][i] = pixD[i];
        _ijk2xyz[i][3] = shift[i];
    }

    create(ndim, imgD, pixD, _ijk2xyz, allocateData);

}

#endif
