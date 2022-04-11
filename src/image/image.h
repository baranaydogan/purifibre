#ifndef SRC_IMAGE_H_
#define SRC_IMAGE_H_

#include "../base/config.h"
#include "niftilib/nifti2/nifti2_io.h"
#include "interpolation.h"


// Features:
// 1. The indexing order of data can be changed by using the "indexData" function.
// e.g. indexData(order), where order is e.g. int order[7] = {3,0,1,2,4,5,6};
// The above example, orders data so that the 4th dimensions run first.
// Note that this does NOT change image dimensions, sub2ind nor ind2sub.
// It also orders the data for faster memory access along a desired dimension.
//
// 2. Interpolation output data type is automatically set based on the caller.
// e.g. Image<bool> img; float a[3]={0.1,0.2,0.3}; float A=img(a); performs the interpolation in float precision.
// i.e. the internally used interpolation coefficients, computations and the return data type are converted to float automatically
// e.g. Image<double> img; double a[3]={0.1,0.2,0.3}; double A=img(a); performs the interpolation in double precision.
// Note that in order to use this functionality, the input point must have the same data type as the output.
// e.g. Image<double> img; float a[3]={0.1,0.2,0.3}; double A=img(a); performs the interpolation in float precision and return a double.

template<typename T>
class Image {

public:
    
    friend class Interpolator<float,T>;
    friend class Interpolator<double,T>;
    friend class Interpolator<long double,T>;
    
    // TODO: Support for complex interpolators
    // friend class Interpolator<std::complex<float>,std::complex<float>>;
    // friend class Interpolator<std::complex<double>,std::complex<double>>;
    // friend class Interpolator<std::complex<long double>,std::complex<long double>>;
    
    Image();
    ~Image();
    Image(std::string _filePath);
    Image(const char* _filePath);
    Image(std::string _filePath, int* _indexOrder);
    Image(const char* _filePath, int* _indexOrder);
    Image(const Image &img);

    std::string   filePath;
    std::string   fileExtension;
    std::string   description;
    SPACEUNIT     spaceUnit;
    TIMEUNIT      timeUnit;
    int           numberOfDimensions;
    int64_t       imgDims[7];
    float         pixDims[7];
    float         xyz2ijk[3][4];
    float         ijk2xyz[3][4];
    DATATYPE      dataType;
    T*            data;

    int64_t       voxCnt; // Number of spatial voxel elements, i.e., voxCnt = imgDims[0]*imgDims[1]*imgDims[2]
    int64_t       valCnt; // Number of data values per voxel elements, i.e., multiplication of all dimensions above 2
    float         smallestPixDim;
    T             outsideVal;    // Value to use for outside the image bound
    int           indexOrder[7];
    
    void          setFilePath(std::string _filePath) {filePath = _filePath; fileExtension = getFileExtension(_filePath);};
    void          setFilePath(const char* _filePath) {setFilePath(std::string(_filePath));};
    bool          readHeader();
    bool          read();
    bool          write(std::string filePath_);
    void          create(int ndims, int64_t* _imgDims, float* _pixDims, float _ijk2xyz[][4], bool allocateData);
    void          createFromTemplate(Image* img, bool allocateData);
    void          createFromTemplate(Image* img, float pixDimScaler, bool allocateData);
    void          createFromBoundingBox(int ndim, std::vector<float> bb, std::vector<float> _pixDims, bool allocateData);    // bb must have values in order, e.g., [x_min x_max y_min y_max z_min z_max t_min t_max ...]
    void          createFromBoundingBox(int ndim, std::vector<float> bb, float _pixDim, bool allocateData);                  // All pixDims are same
    void          createFromBoundingBox(int ndim, std::vector<float> bb, bool allocateData);                                 // Assumes 1 for all pixDims
    void          printInfo();
    
    int           getDimension() {return numberOfDimensions;};
    std::string   getSpaceUnit();
    std::string   getTimeUnit();

    void          allocData()   {parseHeader(); data = new T[voxCnt*valCnt]();   };
    void          deallocData() {if (data!=NULL) { delete[] data; data = NULL; } };
    
    void          setOutsideVal(T _outsideVal) {outsideVal = _outsideVal;}
    T             getOutsideVal() {return outsideVal;}
    
    void          setInterpolationMethod(INTERPMETHOD _interpMethod);
    INTERPMETHOD  getInterpolationMethod() {return interpMethod;}
    
    void          indexData(int* _indexOrder);
    
    int64_t       sub2ind(int64_t i, int64_t j, int64_t k);
    int64_t       sub2ind(int64_t i, int64_t j, int64_t k, int64_t t);
    
    int64_t       sub2ind(int64_t* sub);
    void          ind2sub(int64_t ind, int64_t* sub);
    void          ind2sub(int64_t ind, int64_t& i, int64_t& j, int64_t& k);
    
    template<typename OUT_T>
    void          to_ijk(OUT_T* xyz, OUT_T *ijk);
    
    template<typename OUT_T>
    void          to_xyz(OUT_T*  ijk, OUT_T *xyz);

    template<typename OUT_T>
    void          to_xyz(int64_t ind, OUT_T *xyz);

    bool          isInside(int32_t* ijk);
    bool          isInside(int64_t* ijk);
    bool          isInside(int32_t i, int32_t j, int32_t k);
    bool          isInside(int64_t i, int64_t j, int64_t k);
    
    // Operators
    float       operator()(float *p)                         { return interp3D_f    (this, p);       }
    float       operator()(float *p, int64_t t)              { return interp4D_f_att(this, p, t);    }
    void        operator()(float *p, float* out)             {        interp4D_f_all(this, p, out);  }
    
    double      operator()(double *p)                        { return interp3D_d    (this, p);       }
    double      operator()(double *p, int64_t t)             { return interp4D_d_att(this, p, t);    }
    void        operator()(double *p, double* out)           {        interp4D_d_all(this, p, out);  }
    
    long double operator()(long double *p)                   { return interp3D_ld    (this, p);      }
    long double operator()(long double *p, int64_t t)        { return interp4D_ld_att(this, p, t);   }
    void        operator()(long double *p, long double* out) {        interp4D_ld_all(this, p, out); }

    T           operator()(int64_t i);
    T           operator()(int64_t i, int64_t j, int64_t k);
    void        operator()(int64_t i, int64_t j, int64_t k, T* out);
    
    T*          at(int64_t i, int64_t j, int64_t k);

    // Mgz support
    float         rastkr2ras[4][4];
    
protected:
    
    std::function<float(Image<T>* img, float *p)>                                 interp3D_f;
    std::function<float(Image<T>* img, float *p, int64_t t)>                      interp4D_f_att;
    std::function<void (Image<T>* img, float *p, float* out)>                     interp4D_f_all;
    
    std::function<double(Image<T>* img, double *p)>                               interp3D_d;
    std::function<double(Image<T>* img, double *p, int64_t t)>                    interp4D_d_att;
    std::function<void  (Image<T>* img, double *p, double* out)>                  interp4D_d_all;
    
    std::function<long double(Image<T>* img, long double *p)>                     interp3D_ld;
    std::function<long double(Image<T>* img, long double *p, int64_t t)>          interp4D_ld_att;
    std::function<void       (Image<T>* img, long double *p, long double* out)>   interp4D_ld_all;
    
    INTERPMETHOD    interpMethod;
    
    bool     headerIsRead;
    int64_t  s2i[7];              // multipliers used for sub2ind conversion
    float    dataScaler;
    float    dataOffset;
    DATATYPE inputDataType;       // This is the input data type NOT the data type of this image.
    
private:
    
    void            init();
    void            parseHeader();
    bool            readHeader_nii();
    bool            readHeader_mgz();
    bool            read_nii();
    bool            write_nii(std::string filePath_);
    nifti_1_header* getNiftiHeader();
    
};

// SPEED IS REDUCED IF FUNCTION TEMPLATES ARE DEFINED IN A SEPARATE CPP FILE
template<typename T>
int64_t Image<T>::sub2ind(int64_t i, int64_t j, int64_t k) {
    return i*s2i[0] + j*s2i[1] + k*s2i[2];
}

template<typename T>
int64_t Image<T>::sub2ind(int64_t i, int64_t j, int64_t k, int64_t t) {
    return i*s2i[0] + j*s2i[1] + k*s2i[2] + t*s2i[3];
}

template<typename T>
int64_t Image<T>::sub2ind(int64_t* sub) {
    return sub[0]*s2i[0] + sub[1]*s2i[1] + sub[2]*s2i[2] + sub[3]*s2i[3] + sub[4]*s2i[4] + sub[5]*s2i[5] + sub[6]*s2i[6];
}

template<typename T>
void Image<T>::ind2sub(int64_t ind, int64_t* sub) {
    
    int64_t offset;

    for (int i = 0; i < 7; i++) {
        offset             = ind % imgDims[indexOrder[i]];
        ind               -= offset;
        ind               /= imgDims[indexOrder[i]];
        sub[indexOrder[i]] = offset;
    }
}

template<typename T>
void Image<T>::ind2sub(int64_t ind, int64_t& i, int64_t& j, int64_t& k) {
    i = ind % imgDims[0];
    j = (ind - i)/imgDims[0] % imgDims[1];
    k = ((ind - i)/imgDims[0]-j)/imgDims[1];
}

template<typename T> template<typename OUT_T>
void Image<T>::to_ijk(OUT_T* xyz, OUT_T *ijk) {
    ijk[0]  = xyz2ijk[0][0]*xyz[0] + xyz2ijk[0][1]*xyz[1] + xyz2ijk[0][2]*xyz[2] + xyz2ijk[0][3];
    ijk[1]  = xyz2ijk[1][0]*xyz[0] + xyz2ijk[1][1]*xyz[1] + xyz2ijk[1][2]*xyz[2] + xyz2ijk[1][3];
    ijk[2]  = xyz2ijk[2][0]*xyz[0] + xyz2ijk[2][1]*xyz[1] + xyz2ijk[2][2]*xyz[2] + xyz2ijk[2][3];
}

template<typename T> template<typename OUT_T>
void Image<T>::to_xyz(OUT_T* ijk, OUT_T *xyz) {
    xyz[0]  = ijk2xyz[0][0]*ijk[0] + ijk2xyz[0][1]*ijk[1] + ijk2xyz[0][2]*ijk[2] + ijk2xyz[0][3];
    xyz[1]  = ijk2xyz[1][0]*ijk[0] + ijk2xyz[1][1]*ijk[1] + ijk2xyz[1][2]*ijk[2] + ijk2xyz[1][3];
    xyz[2]  = ijk2xyz[2][0]*ijk[0] + ijk2xyz[2][1]*ijk[1] + ijk2xyz[2][2]*ijk[2] + ijk2xyz[2][3];
}

template<typename T> template<typename OUT_T>
void Image<T>::to_xyz(int64_t ind, OUT_T *xyz) {

    int64_t ijk[3];

    ijk[0]  = ind % imgDims[0];
    ijk[1]  = (ind - ijk[0])/imgDims[0] % imgDims[1];
    ijk[2]  = ((ind - ijk[0])/imgDims[0]-ijk[1])/imgDims[1];

    xyz[0]  = ijk2xyz[0][0]*ijk[0] + ijk2xyz[0][1]*ijk[1] + ijk2xyz[0][2]*ijk[2] + ijk2xyz[0][3];
    xyz[1]  = ijk2xyz[1][0]*ijk[0] + ijk2xyz[1][1]*ijk[1] + ijk2xyz[1][2]*ijk[2] + ijk2xyz[1][3];
    xyz[2]  = ijk2xyz[2][0]*ijk[0] + ijk2xyz[2][1]*ijk[1] + ijk2xyz[2][2]*ijk[2] + ijk2xyz[2][3];
}


template<typename T>
T Image<T>::operator()(int64_t i) {
    return data[i];
}

template<typename T>
T Image<T>::operator()(int64_t i, int64_t j, int64_t k) {
    return ((i<0) || (j<0) || (k<0) || (i>=imgDims[0]) || (j>=imgDims[1]) || (k>=imgDims[2])) ? outsideVal : data[i*s2i[0] + j*s2i[1] + k*s2i[2]];
}

template<typename T>
T* Image<T>::at(int64_t i, int64_t j, int64_t k) {
    return ((i<0) || (j<0) || (k<0) || (i>=imgDims[0]) || (j>=imgDims[1]) || (k>=imgDims[2])) ? NULL : data + i*s2i[0] + j*s2i[1] + k*s2i[2];
}

template<typename T>
void   Image<T>::operator()(int64_t i, int64_t j, int64_t k, T* out) {
    if ((i<0) || (j<0) || (k<0) || (i>=imgDims[0]) || (j>=imgDims[1]) || (k>=imgDims[2]))
        memset(out,outsideVal,valCnt*sizeof(T));
    else
        out = data + i*s2i[0] + j*s2i[1] + k*s2i[2];
}

template<typename T>
void Image<T>::indexData(int* _indexOrder) {
    
    int64_t numel = voxCnt*valCnt;
    
    for (int i=0; i<7; i++) 
        s2i[i] = 1;
    
    for (int i=1; i<7; i++)
        for (int j=0; j<i; j++)
            s2i[_indexOrder[i]] *= imgDims[_indexOrder[j]];
        
    if (data==NULL) {
        for (int i=0; i<7; i++)
            indexOrder[i] = _indexOrder[i];
        return;
    }
    
    int64_t* newIndices = (int64_t*) malloc(numel*sizeof(int64_t));
    
    auto run = [&](MTTASK task) {
        int64_t sub[7];
        ind2sub(task.no,sub);
        newIndices[task.no] = sub[0]*s2i[0] + sub[1]*s2i[1] + sub[2]*s2i[2] + sub[3]*s2i[3] + sub[4]*s2i[4] + sub[5]*s2i[5] + sub[6]*s2i[6];
    };
    MT::MTRUN(numel,MT::maxNumberOfThreads,run);
    
    for (int i=0; i<7; i++)
        indexOrder[i] = _indexOrder[i];
    
    for (int64_t i=0; i<numel; i++) {
        while (newIndices[i] != i) {
            int64_t tmpIndx  = newIndices[newIndices[i]];
            T       tmpData  = data[newIndices[i]];
 
            data[newIndices[i]]       = data[i];
            newIndices[newIndices[i]] = newIndices[i];
 
            newIndices[i] = tmpIndx;
            data[i]       = tmpData;
        }
    }
    
    
    delete[] newIndices;
    
}

template<typename T>
bool Image<T>::isInside(int32_t* ijk) {
    if ( (ijk[0]>-1) && (ijk[1]>-1) && (ijk[2]>-1) && ijk[0]<imgDims[0] && ijk[1]<imgDims[1] && ijk[2]<imgDims[2] )
        return true;
    else 
        return false;
}

template<typename T>
bool Image<T>::isInside(int64_t* ijk) {
    if ( (ijk[0]>-1) && (ijk[1]>-1) && (ijk[2]>-1) && ijk[0]<imgDims[0] && ijk[1]<imgDims[1] && ijk[2]<imgDims[2] )
        return true;
    else 
        return false;
}

template<typename T>
bool Image<T>::isInside(int32_t i, int32_t j, int32_t k) {
    if ( (i>-1) && (j>-1) && (k>-1) && i<imgDims[0] && j<imgDims[1] && k<imgDims[2] )
        return true;
    else 
        return false;
}

template<typename T>
bool Image<T>::isInside(int64_t i, int64_t j, int64_t k) {
    if ( (i>-1) && (j>-1) && (k>-1) && i<imgDims[0] && j<imgDims[1] && k<imgDims[2] )
        return true;
    else 
        return false;
}

template<typename T>
void   Image<T>::setInterpolationMethod(INTERPMETHOD _interpMethod) {
    
    interpMethod = _interpMethod; 
    
    switch (interpMethod) {
        case NEAREST:
            interp3D_f       = &Interpolator<float,      T>::interp_nearest_3D;
            interp4D_f_att   = &Interpolator<float,      T>::interp_nearest_4D_att;
            interp4D_f_all   = &Interpolator<float,      T>::interp_nearest_4D_all;
            interp3D_d       = &Interpolator<double,     T>::interp_nearest_3D;
            interp4D_d_att   = &Interpolator<double,     T>::interp_nearest_4D_att;
            interp4D_d_all   = &Interpolator<double,     T>::interp_nearest_4D_all;
            interp3D_ld      = &Interpolator<long double,T>::interp_nearest_3D;
            interp4D_ld_att  = &Interpolator<long double,T>::interp_nearest_4D_att;
            interp4D_ld_all  = &Interpolator<long double,T>::interp_nearest_4D_all;
            break;
        case LINEAR:
            interp3D_f       = &Interpolator<float,      T>::interp_linear_3D;
            interp4D_f_att   = &Interpolator<float,      T>::interp_linear_4D_att;
            interp4D_f_all   = &Interpolator<float,      T>::interp_linear_4D_all;
            interp3D_d       = &Interpolator<double,     T>::interp_linear_3D;
            interp4D_d_att   = &Interpolator<double,     T>::interp_linear_4D_att;
            interp4D_d_all   = &Interpolator<double,     T>::interp_linear_4D_all;
            interp3D_ld      = &Interpolator<long double,T>::interp_linear_3D;
            interp4D_ld_att  = &Interpolator<long double,T>::interp_linear_4D_att;
            interp4D_ld_all  = &Interpolator<long double,T>::interp_linear_4D_all;
            break;
        case CUBIC:
            interp3D_f       = &Interpolator<float,      T>::interp_cubic_3D;
            interp4D_f_att   = &Interpolator<float,      T>::interp_cubic_4D_att;
            interp4D_f_all   = &Interpolator<float,      T>::interp_cubic_4D_all;
            interp3D_d       = &Interpolator<double,     T>::interp_cubic_3D;
            interp4D_d_att   = &Interpolator<double,     T>::interp_cubic_4D_att;
            interp4D_d_all   = &Interpolator<double,     T>::interp_cubic_4D_all;
            interp3D_ld      = &Interpolator<long double,T>::interp_cubic_3D;
            interp4D_ld_att  = &Interpolator<long double,T>::interp_cubic_4D_att;
            interp4D_ld_all  = &Interpolator<long double,T>::interp_cubic_4D_all;
            break;
        default:
            interp3D_f       = &Interpolator<float,      T>::interp_linear_3D;
            interp4D_f_att   = &Interpolator<float,      T>::interp_linear_4D_att;
            interp4D_f_all   = &Interpolator<float,      T>::interp_linear_4D_all;
            interp3D_d       = &Interpolator<double,     T>::interp_linear_3D;
            interp4D_d_att   = &Interpolator<double,     T>::interp_linear_4D_att;
            interp4D_d_all   = &Interpolator<double,     T>::interp_linear_4D_all;
            interp3D_ld      = &Interpolator<long double,T>::interp_linear_3D;
            interp4D_ld_att  = &Interpolator<long double,T>::interp_linear_4D_att;
            interp4D_ld_all  = &Interpolator<long double,T>::interp_linear_4D_all;
            break;
    }
    
    
}



#endif
