#ifndef SRC_IMAGE_READER_CPP_
#define SRC_IMAGE_READER_CPP_

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
bool Image<T>::read() {
    
    if (headerIsRead==false) {
        if (readHeader()==false)
            return false;
    }
    
    if ((fileExtension!="nii.gz") && (fileExtension!="nii")) {
        std::cout << "Can't read image data with this extension yet: " << fileExtension << std::endl << std::flush;
        return false;
    } else {
        return read_nii();
    }
    
    return true;
    
}

template <typename OUT_T,typename INP_T>
void convert(OUT_T* out, void* inp, int64_t* imgDims, float dataScaler, float dataOffset)
{                

    int64_t numel = 1;
    for (auto i=0; i<7; i++) numel *= imgDims[i];
    
    MT::MTRUN(numel,MT::maxNumberOfThreads,[&](MTTASK task) {out[task.no] = *((INP_T*)inp+task.no);});
    
    bool scaleData = ( (dataScaler==0) || ((dataScaler==1) && (dataOffset==0)) ) ? false : true;
    if (scaleData) {
        MT::MTRUN(numel,MT::maxNumberOfThreads,[&](MTTASK task) {out[task.no] = dataScaler*out[task.no] + dataOffset;});
    }
    
}

template <typename OUT_T,typename INP_T>
void convert(OUT_T* out, void* inp, int64_t* imgDims, int* outIndexOrder, float dataScaler, float dataOffset)
{            
    int64_t outS2i[7];
    
    for (int i=0; i<7; i++) 
        outS2i[i] = 1;
    
    for (int i=1; i<7; i++)
        for (int j=0; j<i; j++)
            outS2i[outIndexOrder[i]] *= imgDims[outIndexOrder[j]];
    
    int64_t numel = 1;
    for (auto i=0; i<7; i++) numel *= imgDims[i];
    
    auto indexData = [&](MTTASK task) {
        
        int64_t sub[7];
        int64_t offset;
        int64_t index = task.no;

        for (int i = 0; i < 7; i++) {
            offset     = index % imgDims[i];
            index     -= offset;
            index     /= imgDims[i];
            sub[i]     = offset;
        }

        out[sub[0]*outS2i[0] + sub[1]*outS2i[1] + sub[2]*outS2i[2] + sub[3]*outS2i[3] + sub[4]*outS2i[4] + sub[5]*outS2i[5] + sub[6]*outS2i[6]] = *((INP_T*)inp+task.no);
        
    };
    MT::MTRUN(numel,MT::maxNumberOfThreads,indexData);
    
    
    bool scaleData = ( (dataScaler==0) || ((dataScaler==1) && (dataOffset==0)) ) ? false : true;
    if (scaleData) {
        MT::MTRUN(numel,MT::maxNumberOfThreads,[&](MTTASK task) {out[task.no] = dataScaler*out[task.no] + dataOffset;});
    }
    
}

template<typename T>
bool Image<T>::read_nii() {
    
    nifti_image* nim = nifti_image_read(filePath.c_str(),0);
    
    if (nifti_image_load(nim)==-1) {
        std::cout<<"Cannot read nifti image: " << filePath << std::endl;
        return false;
    }
    
    int64_t numel = voxCnt*valCnt;
    
    data = (T*) malloc(numel*sizeof(T));
    
    bool reIdx = false;
    for (int i=1; i<7; i++)
        if (indexOrder[i] != i)
            reIdx = true;
    
    switch (inputDataType) {
        
        case UINT8:         reIdx ? convert<T,uint8_t>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,uint8_t>    (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case INT8:          reIdx ? convert<T,int8_t>     (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,int8_t>     (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case UINT16:        reIdx ? convert<T,uint16_t>   (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,uint16_t>   (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case INT16:         reIdx ? convert<T,int16_t>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,int16_t>    (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case UINT32:        reIdx ? convert<T,uint32_t>   (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,uint32_t>   (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case INT32:         reIdx ? convert<T,int32_t>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,int32_t>    (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case UINT64:        reIdx ? convert<T,uint64_t>   (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,uint64_t>   (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case INT64:         reIdx ? convert<T,int64_t>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,int64_t>    (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case FLOAT32:       reIdx ? convert<T,float>      (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,float>      (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case FLOAT64:       reIdx ? convert<T,double>     (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,double>     (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        case FLOAT128:      reIdx ? convert<T,long double>(data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,long double>(data,nim->data,imgDims,dataScaler,dataOffset);    break;
        
        // TODO: Implement converters for complex data types
        // case COMPLEX64:     reIdx ? convert<T,std::complex<double>)>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,std::complex<double>>    (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        // case COMPLEX128:    reIdx ? convert<T,std::complex<long double>)>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,std::complex<long double>>    (data,nim->data,imgDims,dataScaler,dataOffset);    break;
        // case COMPLEX256:    reIdx ? convert<T,std::complex<long long double>)>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset) : convert<T,std::complex<long long double>>    (data,nim->data,imgDims,dataScaler,dataOffset);    break;
            
        default:   
            std::cout<<"Can't read nifti file. Unknown datatype" << std::endl; 
            break;
    }
    
    nifti_image_free(nim);
    
    return true;
}

#endif
