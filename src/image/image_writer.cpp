#ifndef SRC_IMAGE_WRITER_CPP_
#define SRC_IMAGE_WRITER_CPP_

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
bool Image<T>::write(std::string filePath_) {
    
    std::string ext = getFileExtension(filePath_);
    
    if ((ext=="nii.gz") || (ext=="nii"))
        return write_nii(filePath_);
    
    return true;
    
}

template<typename T>
nifti_1_header* Image<T>::getNiftiHeader() {
    
    nifti_1_header* header; 
    
    const int64_t dims[8] = {numberOfDimensions,imgDims[0],imgDims[1],imgDims[2],imgDims[3],imgDims[4],imgDims[5],imgDims[6]};
    
    switch (dataType) {
    case BOOL:          header = nifti_make_new_n1_header(dims,DT_UINT8    ); break;
    case UINT8:         header = nifti_make_new_n1_header(dims,DT_UINT8    ); break;
    case INT8:          header = nifti_make_new_n1_header(dims,DT_INT8     ); break;
    case UINT16:        header = nifti_make_new_n1_header(dims,DT_UINT16   ); break;
    case INT16:         header = nifti_make_new_n1_header(dims,DT_INT16    ); break;
    case UINT32:        header = nifti_make_new_n1_header(dims,DT_UINT32   ); break;
    case INT32:         header = nifti_make_new_n1_header(dims,DT_INT32    ); break;
    case UINT64:        header = nifti_make_new_n1_header(dims,DT_UINT64   ); break;
    case INT64:         header = nifti_make_new_n1_header(dims,DT_INT64    ); break;
    case FLOAT32:       header = nifti_make_new_n1_header(dims,DT_FLOAT32  ); break;
    case FLOAT64:       header = nifti_make_new_n1_header(dims,DT_FLOAT64  ); break;
    case FLOAT128:      header = nifti_make_new_n1_header(dims,DT_FLOAT128 ); break;
    default:   msg_error("Unknown output datatype."); return NULL; break;
    }
    
    header->scl_slope       = 1;
    header->scl_inter       = 0;

    
    header->qform_code      = 1;
    header->sform_code      = 1;
    
    
    nifti_dmat44 R;
    for (int i=0; i<3; i++)
        for (int j=0; j<4; j++) {
            R.m[i][j] = ijk2xyz[i][j];
        }
    R.m[3][0] = R.m[3][1] = R.m[3][2] = 0;
    R.m[3][3] = 1;
    
    
    double qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac; // not used
    nifti_dmat44_to_quatern(R,&qb,&qc,&qd,&qx,&qy,&qz,&dx,&dy,&dz,&qfac);
    
    header->quatern_b = qb;
    header->quatern_c = qc;
    header->quatern_d = qd;
    header->qoffset_x = qx;
    header->qoffset_y = qy;
    header->qoffset_z = qz;
    header->pixdim[0] = qfac;
        
    header->dim[0] = numberOfDimensions;
    for (int i=0; i<numberOfDimensions; i++) {
        header->pixdim[i+1] = pixDims[i];
        header->dim[i+1]    = imgDims[i];
    }
    
    header->srow_x[0] = R.m[0][0];
    header->srow_x[1] = R.m[0][1];
    header->srow_x[2] = R.m[0][2];
    header->srow_x[3] = R.m[0][3];
    header->srow_y[0] = R.m[1][0];
    header->srow_y[1] = R.m[1][1];
    header->srow_y[2] = R.m[1][2];
    header->srow_y[3] = R.m[1][3];
    header->srow_z[0] = R.m[2][0];
    header->srow_z[1] = R.m[2][1];
    header->srow_z[2] = R.m[2][2];
    header->srow_z[3] = R.m[2][3];
    
    switch (spaceUnit) {
    case METER:     header->xyzt_units = NIFTI_UNITS_METER;      break;
    case MM:        header->xyzt_units = NIFTI_UNITS_MM;         break;
    case MICRON:    header->xyzt_units = NIFTI_UNITS_MICRON;     break;
    default:        header->xyzt_units = NIFTI_UNITS_UNKNOWN;    break;
    }
    
    for (size_t i=0; i<SGNTR.length() && i<80; i++) {
        header->descrip[i] = SGNTR[i];
    }
    
    return header;
    
}


template<typename T>
void* resetIndexingAndCopyData(T* inp, int64_t* imgDims, int* indexOrder) {
 
    int64_t s2i[7]; 
    
    int64_t numel = 1;
    
    for (int i=0; i<7; i++) { 
        s2i[i] = 1;
        numel *= imgDims[i];
    }
    
    for (int i=1; i<7; i++)
        for (int j=0; j<i; j++)
            s2i[i] *= imgDims[j];
    
    void* out = (T*) malloc(numel*sizeof(T));
        
    auto run = [&](MTTASK task) {
        int64_t sub[7];
        int64_t offset;
        int64_t ind = task.no;

        for (int i = 0; i < 7; i++) {
            offset             = ind % imgDims[indexOrder[i]];
            ind               -= offset;
            ind               /= imgDims[indexOrder[i]];
            sub[indexOrder[i]] = offset;
        }
        
        *((T*)(out)+sub[0]*s2i[0] + sub[1]*s2i[1] + sub[2]*s2i[2] + sub[3]*s2i[3] + sub[4]*s2i[4] + sub[5]*s2i[5] + sub[6]*s2i[6]) = inp[task.no];
        
    };
    MT::MTRUN(numel,MT::maxNumberOfThreads,run);
    
    return out;

}

template<typename T>
bool Image<T>::write_nii(std::string filePath_) {
    
    nifti_1_header* header = getNiftiHeader();
    nifti_image* nim       = nifti_convert_n1hdr2nim(*header,filePath_.c_str());
    
    bool resetIndexing = false;
    
    for (int i=0; i<7; i++)
        if (indexOrder[i]!=i)
            resetIndexing = true;
    
    if (resetIndexing)
        nim->data = resetIndexingAndCopyData(data,imgDims,indexOrder);
    else
        nim->data = data;
    
    nifti_image_write(nim);
    
    if (resetIndexing) {
        nifti_image_free(nim);
    } else {
        if( nim->fname != NULL ) free(nim->fname) ;
        if( nim->iname != NULL ) free(nim->iname) ;
        (void)nifti_free_extensions( nim ) ;
        free(nim);
    }
    
    free(header);
    return true;
    
}

#endif
