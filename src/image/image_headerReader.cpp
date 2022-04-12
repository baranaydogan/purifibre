#ifndef SRC_IMAGE_HEADERREADER_CPP_
#define SRC_IMAGE_HEADERREADER_CPP_

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
bool Image<T>::readHeader() {

    if (headerIsRead) {

        return true;

    } else {

        if ((fileExtension=="nii.gz") || (fileExtension=="nii"))
            return readHeader_nii();

        if (fileExtension=="mgz")
            return readHeader_mgz();

        std::cout << "Unknown file extension: " << fileExtension << std::endl << std::flush;
        return false;

    }

}


template<typename T>
void Image<T>::parseHeader() {

    voxCnt   = imgDims[0]*imgDims[1]*imgDims[2];
    valCnt   = imgDims[3]*imgDims[4]*imgDims[5]*imgDims[6];

    if (valCnt<1) valCnt = 1;

    smallestPixDim = pixDims[0];
    smallestPixDim = (smallestPixDim < pixDims[1]) ? smallestPixDim : pixDims[1];
    smallestPixDim = (smallestPixDim < pixDims[2]) ? smallestPixDim : pixDims[2];

    // indexOrder and s2i are used for indexing the data in any desired order and for sub2ind and ind2sub conversion
    // These do not change the image dimensions.

    for (int i=0; i<7; i++)
        s2i[i] = 1;

    for (int i=1; i<7; i++)
        for (int j=0; j<i; j++)
            s2i[indexOrder[i]] *= imgDims[indexOrder[j]];

    headerIsRead   = true;

}

template<typename T>
bool Image<T>::readHeader_nii() {

    nifti_image* nim = nifti_image_read(filePath.c_str(),0);
    if (nim==NULL) exit(EXIT_FAILURE);

    description = std::string(nim->descrip);

    switch (nim->xyz_units) {
    case NIFTI_UNITS_METER:     spaceUnit = METER;              break;
    case NIFTI_UNITS_MM:        spaceUnit = MM;                 break;
    case NIFTI_UNITS_MICRON:    spaceUnit = MICRON;             break;
    default:                    spaceUnit = UNKNOWNSPACEUNIT;   break;
    }

    switch (nim->time_units) {
    case NIFTI_UNITS_SEC:       timeUnit = SEC;             break;
    case NIFTI_UNITS_MSEC:      timeUnit = MSEC;            break;
    case NIFTI_UNITS_USEC:      timeUnit = USEC;            break;
    case NIFTI_UNITS_HZ:        timeUnit = HZ;              break;
    case NIFTI_UNITS_PPM:       timeUnit = PPM;             break;
    case NIFTI_UNITS_RADS:      timeUnit = RADS;            break;
    default:                    timeUnit = UNKNOWNTIMEUNIT; break;
    }

    dataScaler = nim->scl_slope;
    dataOffset = nim->scl_inter;

    switch (nim->datatype) {

    case DT_UINT8:         inputDataType=UINT8_DT;      break;
    case DT_INT8:          inputDataType=INT8_DT;       break;
    case DT_UINT16:        inputDataType=UINT16_DT;     break;
    case DT_INT16:         inputDataType=INT16_DT;      break;
    case DT_UINT32:        inputDataType=UINT32_DT;     break;
    case DT_INT32:         inputDataType=INT32_DT;      break;
    case DT_UINT64:        inputDataType=UINT64_DT;     break;
    case DT_INT64:         inputDataType=INT64_DT;      break;
    case DT_FLOAT32:       inputDataType=FLOAT32_DT;    break;
    case DT_FLOAT64:       inputDataType=FLOAT64_DT;    break;
    case DT_FLOAT128:      inputDataType=FLOAT128_DT;   break;

    case DT_COMPLEX64:
        inputDataType=COMPLEX64_DT;
        std::cout<<"Nifti datatype: complex64 "       << std::endl; std::cout << "is not an accepted datatype" << std::endl;
        break;

    case DT_COMPLEX128:
        inputDataType=COMPLEX128_DT;
        std::cout<<"Nifti datatype: complex128 "      << std::endl; std::cout << "is not an accepted datatype" << std::endl;
        break;

    case DT_COMPLEX256:
        inputDataType=COMPLEX256_DT;
        std::cout<<"Nifti datatype: complex256 "      << std::endl; std::cout << "is not an accepted datatype" << std::endl;
        break;

    case DT_BINARY:
        inputDataType=UNKNOWN_DT;
        std::cout<<"Nifti datatype: binary "          << std::endl; std::cout << "is not an accepted datatype" << std::endl;
        break;

    case DT_RGB:
        inputDataType=UNKNOWN_DT;
        std::cout<<"Nifti datatype: rgb24 "           << std::endl; std::cout << "is not an accepted datatype" << std::endl;
        break;

    case DT_RGBA32:
        inputDataType=UNKNOWN_DT;
        std::cout<<"Nifti datatype: rgba32 "          << std::endl; std::cout << "is not an accepted datatype" << std::endl;
        break;

    case DT_ALL:
        inputDataType=UNKNOWN_DT;
        std::cout<<"Nifti datatype: all  "            << std::endl; std::cout << "is not an accepted datatype" << std::endl;
        break;

    default:
        inputDataType=UNKNOWN_DT; 
        std::cout<<"Nifti datatype: unknown or not applicable "        << std::endl; std::cout << "is not an accepted datatype" << std::endl;
        break;
    }

    // Get dims and pixDims
    numberOfDimensions = nim->dim[0];

    if (numberOfDimensions>0) {
        for (int i=0; i<7; i++) {
            imgDims[i] = nim->dim[i+1];
            pixDims[i] = (imgDims[i]==0) ? 1 : nim->pixdim[i+1];
            imgDims[i] = (imgDims[i]==0) ? 1 : imgDims[i];
        }
    }

    // Use sform if possible otherwise use qform
    if (nim->sform_code>0) {
        for (int i=0; i<3; i++)
            for (int j=0; j<4; j++) {
                xyz2ijk[i][j] = nim->sto_ijk.m[i][j];
                ijk2xyz[i][j] = nim->sto_xyz.m[i][j];
            }
    }
    else {
        for (int i=0; i<3; i++)
            for (int j=0; j<4; j++) {
                xyz2ijk[i][j] = nim->qto_ijk.m[i][j];
                ijk2xyz[i][j] = nim->qto_xyz.m[i][j];
            }
    }

    nifti_image_free(nim);

    parseHeader();

    return true;

}

template<typename T>
bool Image<T>::readHeader_mgz() {

    size_t posExt     = filePath.find_last_of(".");
    size_t posPath    = filePath.find_last_of("/");

    std::string path      = filePath.substr(0,posPath);
    std::string base      = filePath.substr(posPath+1, posExt-posPath-1);
    std::string extension = filePath.substr(posExt+1);

    if (extension == "mgz") {

        char decomFname[L_tmpnam];

        // TODO: review if the following solution is working as intended
        // if (mkstemp(decomFname)) {
        if (std::tmpnam(decomFname)) {
            std::string unzipCommand = "gunzip -c " + filePath + " > " + std::string(decomFname);
            system(unzipCommand.c_str());
        } else  {
            std::cout << "Can't create tmp file " << extension << std::endl;
            return false;
        }

        FILE *inp;
        inp = fopen(decomFname,"r");

        // int version;
        // int type;
        // int dof;
        short goodRASFlag;
        float xr,xa,xs;
        float yr,ya,ys;
        float zr,za,zs;
        float cr,ca,cs;

        int   tmpi;
        short tmps;
        float tmpf;
        std::fread(&tmpi,sizeof(int),  1,inp); swapByteOrder(tmpi);   // version     = tmpi;
        std::fread(&tmpi,sizeof(int),  1,inp); swapByteOrder(tmpi);   imgDims[0]     = tmpi;
        std::fread(&tmpi,sizeof(int),  1,inp); swapByteOrder(tmpi);   imgDims[1]     = tmpi;
        std::fread(&tmpi,sizeof(int),  1,inp); swapByteOrder(tmpi);   imgDims[2]     = tmpi;
        std::fread(&tmpi,sizeof(int),  1,inp); swapByteOrder(tmpi);   imgDims[3]     = tmpi;
        numberOfDimensions = (imgDims[3]>1) ? 4 : 3;

        std::fread(&tmpi,sizeof(int),  1,inp); swapByteOrder(tmpi);   // type        = tmpi;
        std::fread(&tmpi,sizeof(int),  1,inp); swapByteOrder(tmpi);   // dof         = tmpi;
        std::fread(&tmps,sizeof(short),1,inp); swapByteOrder(tmps); goodRASFlag = tmps;

        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); pixDims[0]  = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); pixDims[1]  = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); pixDims[2]  = tmpf;

        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); xr = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); xa = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); xs = tmpf;

        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); yr = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); ya = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); ys = tmpf;

        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); zr = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); za = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); zs = tmpf;

        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); cr = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); ca = tmpf;
        std::fread(&tmpf,sizeof(float),1,inp); swapByteOrder(tmpf); cs = tmpf;

        if (goodRASFlag!=1) {
            std::cout << "File is not good for RAS mm conversion." << std::endl;
            return false;
        }



        // std::cout << "Version    : " << version     << std::endl;
        // std::cout << "Width      : " << dims[0]     << std::endl;
        // std::cout << "Height     : " << dims[1]     << std::endl;
        // std::cout << "Depth      : " << dims[2]     << std::endl;
        // std::cout << "Volumes    : " << dims[3]     << std::endl;

        // std::cout << "type       : " << type        << std::endl;
        // std::cout << "dof        : " << dof         << std::endl;
        // std::cout << "goodRASFlag: " << goodRASFlag << std::endl;
        // std::cout << "xdim       : " << pixDims[0]  << std::endl;
        // std::cout << "ydim       : " << pixDims[1]  << std::endl;
        // std::cout << "zdim       : " << pixDims[2]  << std::endl;

        // std::cout << "xr         : " << xr          << std::endl;
        // std::cout << "xa         : " << xa          << std::endl;
        // std::cout << "xs         : " << xs          << std::endl;
        // std::cout << "yr         : " << yr          << std::endl;
        // std::cout << "ya         : " << ya          << std::endl;
        // std::cout << "ys         : " << ys          << std::endl;
        // std::cout << "zr         : " << zr          << std::endl;
        // std::cout << "za         : " << za          << std::endl;
        // std::cout << "zs         : " << zs          << std::endl;
        // std::cout << "cr         : " << cr          << std::endl;
        // std::cout << "ca         : " << ca          << std::endl;
        // std::cout << "cs         : " << cs          << std::endl;

        // Not used
        // float xstart = -dims[0]/2.0*pixDims[0];
        // float xend   =  xstart;
        // float ystart = -dims[1]/2.0*pixDims[1];
        // float yend   =  ystart;
        // float zstart = -dims[2]/2.0*pixDims[2];
        // float zend   =  zstart;


        // Choose between sform or qform
        float ci    = imgDims[0]/2.0;
        float cj    = imgDims[1]/2.0;
        float ck    = imgDims[2]/2.0;

        mat44 vox2ras;
        vox2ras.m[0][0] = pixDims[0]*xr;
        vox2ras.m[0][1] = pixDims[1]*yr;
        vox2ras.m[0][2] = pixDims[2]*zr;
        vox2ras.m[0][3] = cr - (vox2ras.m[0][0]*ci + vox2ras.m[0][1]*cj + vox2ras.m[0][2]*ck);
        vox2ras.m[1][0] = pixDims[0]*xa;
        vox2ras.m[1][1] = pixDims[1]*ya;
        vox2ras.m[1][2] = pixDims[2]*za;
        vox2ras.m[1][3] = ca - (vox2ras.m[1][0]*ci + vox2ras.m[1][1]*cj + vox2ras.m[1][2]*ck);
        vox2ras.m[2][0] = pixDims[0]*xs;
        vox2ras.m[2][1] = pixDims[1]*ys;
        vox2ras.m[2][2] = pixDims[2]*zs;
        vox2ras.m[2][3] = cs - (vox2ras.m[2][0]*ci + vox2ras.m[2][1]*cj + vox2ras.m[2][2]*ck);
        vox2ras.m[3][0] = 0;
        vox2ras.m[3][1] = 0;
        vox2ras.m[3][2] = 0;
        vox2ras.m[3][3] = 0;

        mat44 vox2rastkr;
        vox2rastkr.m[0][0] = pixDims[0]*xr;
        vox2rastkr.m[0][1] = pixDims[1]*yr;
        vox2rastkr.m[0][2] = pixDims[2]*zr;
        vox2rastkr.m[0][3] = -(vox2ras.m[0][0]*ci + vox2ras.m[0][1]*cj + vox2ras.m[0][2]*ck);
        vox2rastkr.m[1][0] = pixDims[0]*xa;
        vox2rastkr.m[1][1] = pixDims[1]*ya;
        vox2rastkr.m[1][2] = pixDims[2]*za;
        vox2rastkr.m[1][3] = -(vox2ras.m[1][0]*ci + vox2ras.m[1][1]*cj + vox2ras.m[1][2]*ck);
        vox2rastkr.m[2][0] = pixDims[0]*xs;
        vox2rastkr.m[2][1] = pixDims[1]*ys;
        vox2rastkr.m[2][2] = pixDims[2]*zs;
        vox2rastkr.m[2][3] = -(vox2ras.m[2][0]*ci + vox2ras.m[2][1]*cj + vox2ras.m[2][2]*ck);
        vox2rastkr.m[3][0] = 0;
        vox2rastkr.m[3][1] = 0;
        vox2rastkr.m[3][2] = 0;
        vox2rastkr.m[3][3] = 0;

        mat44 rastkr2vox;
        rastkr2vox = nifti_mat44_inverse(vox2rastkr);


        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                rastkr2ras[i][j] =
                vox2ras.m[i][0]*rastkr2vox.m[0][j] +
                vox2ras.m[i][1]*rastkr2vox.m[1][j] +
                vox2ras.m[i][2]*rastkr2vox.m[2][j] +
                vox2ras.m[i][3]*rastkr2vox.m[3][j];
            }
        }

        fclose(inp);

        std::string deleteTmpFileCommand = "rm " + std::string(decomFname);
        system(deleteTmpFileCommand.c_str());


    } else {
        std::cout << "Can't read file with extension " << extension << std::endl;
        return false;
    }

    parseHeader();

    // TODO:
    // Make xyz2ijk and ijk2xyz
    // Make a nim here based on the above information
    // const int dims[8] = {3,10,8,20,1,1,1,1};
    // nifti_make_new_nim(dims, DT_UINT32, true);
    // nim->datatype = DT_UINT32;


    return true;

}

#endif
