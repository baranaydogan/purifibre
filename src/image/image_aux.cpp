#ifndef SRC_IMAGE_AUX_CPP_
#define SRC_IMAGE_AUX_CPP_

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
std::string Image<T>::getSpaceUnit() {
    std::string unit;
    switch (spaceUnit) {
    case METER:         unit = "meter";        break;
    case MM:            unit = "mm";           break;
    case MICRON:        unit = "micron";       break;
    default:            unit = "unknown unit"; break;
    }
    return unit;
}

template<typename T>
void Image<T>::printInfo() {

    std::cout << std::endl;
    std::cout << "File name:              "   << filePath << std::endl;
    std::cout << "File type:              "   << fileExtension << std::endl;
    std::cout << "Data description:       "   << description << std::endl;
    std::cout << "Number of dimensions:   "   << numberOfDimensions << std::endl;
    std::cout << "Dimensions:             [";
    for (int i=0; i<numberOfDimensions; i++) {
        std::cout << imgDims[i];
        if (i!=(numberOfDimensions-1))
            std::cout << " x ";
    }
    std::cout << "]" << std::endl;
    std::cout << "Pixdim:                 [" << pixDims[0] << " x " << pixDims[1] << " x " << pixDims[2] << "]" << std::endl;
    std::cout << "Unit of pixdim:         " << getSpaceUnit() << std::endl;

    std::cout << "Datatype:               ";

    switch (inputDataType) {

        case UINT8_DT:         std::cout << "UINT8";      break;
        case INT8_DT:          std::cout << "INT8";       break;
        case UINT16_DT:        std::cout << "UINT16";     break;
        case INT16_DT:         std::cout << "INT16";      break;
        case UINT32_DT:        std::cout << "UINT32";     break;
        case INT32_DT:         std::cout << "INT32";      break;
        case UINT64_DT:        std::cout << "UINT64";     break;
        case INT64_DT:         std::cout << "INT64";      break;
        case FLOAT32_DT:       std::cout << "FLOAT32";    break;
        case FLOAT64_DT:       std::cout << "FLOAT64";    break;
        case FLOAT128_DT:      std::cout << "FLOAT128";   break;

        default:
            std::cout<<"Unknown data type";
        break;
    }

    std::cout << std::endl;

    std::cout << "Datatype (internal use):" << TYPENAMES[typeid(T)] << std::endl;

    std::cout << std::setprecision(4) << "ijk2xyz:                " << std::endl;
    std::cout << std::setprecision(4) << "                        " << ijk2xyz[0][0] << " " << ijk2xyz[0][1] << " " << ijk2xyz[0][2] << " " << ijk2xyz[0][3] << std::endl;
    std::cout << std::setprecision(4) << "                        " << ijk2xyz[1][0] << " " << ijk2xyz[1][1] << " " << ijk2xyz[1][2] << " " << ijk2xyz[1][3] << std::endl;
    std::cout << std::setprecision(4) << "                        " << ijk2xyz[2][0] << " " << ijk2xyz[2][1] << " " << ijk2xyz[2][2] << " " << ijk2xyz[2][3] << std::endl;

    std::cout << std::setprecision(4) << "xyz2ijk:                " << std::endl;
    std::cout << std::setprecision(4) << "                        " << xyz2ijk[0][0] << " " << xyz2ijk[0][1] << " " << xyz2ijk[0][2] << " " << xyz2ijk[0][3] << std::endl;
    std::cout << std::setprecision(4) << "                        " << xyz2ijk[1][0] << " " << xyz2ijk[1][1] << " " << xyz2ijk[1][2] << " " << xyz2ijk[1][3] << std::endl;
    std::cout << std::setprecision(4) << "                        " << xyz2ijk[2][0] << " " << xyz2ijk[2][1] << " " << xyz2ijk[2][2] << " " << xyz2ijk[2][3] << std::endl;

}

#endif
