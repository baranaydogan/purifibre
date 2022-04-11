#include "tractogram2imageMapper.h"
#include <set>

// Explicit instantiations
template class Tractogram2ImageMapper<bool>;
template class Tractogram2ImageMapper<uint8_t>;
template class Tractogram2ImageMapper<int8_t>;
template class Tractogram2ImageMapper<uint16_t>;
template class Tractogram2ImageMapper<int16_t>;
template class Tractogram2ImageMapper<uint32_t>;
template class Tractogram2ImageMapper<int32_t>;
template class Tractogram2ImageMapper<uint64_t>;
template class Tractogram2ImageMapper<int64_t>;
template class Tractogram2ImageMapper<float>;
template class Tractogram2ImageMapper<double>;
template class Tractogram2ImageMapper<long double>;

// TODO: These are not yet supported because the Image<T> is missing their definitions
// template class Tractogram2ImageMapper<std::complex<float>>;
// template class Tractogram2ImageMapper<std::complex<double>>;
// template class Tractogram2ImageMapper<std::complex<long double>>;

template<typename T>
Tractogram2ImageMapper<T>::Tractogram2ImageMapper(TractogramReader* _tractogram, Image<T>* _img) {

    // Initialize tractogram and make copies for multithreader
    tractogram = new TractogramReader[MT::maxNumberOfThreads]();
    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        tractogram[t].copyFrom(*_tractogram);
    }

    mask       = NULL;
    mapOnce    = false;
    img        = _img;

    img->readHeader();

    grid.resize(MT::maxNumberOfThreads);
    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        grid[t].resize(img->imgDims[0]);
        for (int i = 0; i < img->imgDims[0]; i++) {
            grid[t][i].resize(img->imgDims[1]);
            for (int j = 0; j < img->imgDims[1]; j++) {
                grid[t][i][j].resize(img->imgDims[2]);
            }
        }
    }

    gridAllocated = true;
    smoothing = std::make_tuple(0,0);

}

template<typename T>
Tractogram2ImageMapper<T>::Tractogram2ImageMapper(TractogramReader* _tractogram, Image<T>* _img, bool allocateGrid) {

    // Initialize tractogram and make copies for multithreader
    tractogram = new TractogramReader[MT::maxNumberOfThreads]();
    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        tractogram[t].copyFrom(*_tractogram);
    }


    mask       = NULL;
    mapOnce    = false;
    img        = _img;

    img->readHeader();

    if (allocateGrid) {
        grid.resize(MT::maxNumberOfThreads);
        for (int t = 0; t < MT::maxNumberOfThreads; t++) {
            grid[t].resize(img->imgDims[0]);
            for (int i = 0; i < img->imgDims[0]; i++) {
                grid[t][i].resize(img->imgDims[1]);
                for (int j = 0; j < img->imgDims[1]; j++) {
                    grid[t][i][j].resize(img->imgDims[2],NULL);
                }
            }
        }
    }

    gridAllocated = allocateGrid;
    smoothing = std::make_tuple(0,0);

}

template<typename T>
Tractogram2ImageMapper<T>::~Tractogram2ImageMapper() {

    if (mask!=NULL) {
        for (int i = 0; i < img->imgDims[0]; i++) {
            for (int j = 0; j < img->imgDims[1]; j++) {
                delete[] mask[i][j];
            }
            delete[] mask[i];
        }
        delete[] mask;
    }

    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        tractogram[t].destroyCopy();
    }
    delete[] tractogram;

}

template<typename T>
bool Tractogram2ImageMapper<T>::setMask(Image<int>* maskImg) {
        
    maskImg->readHeader();

    for (int i=0; i<3; i++) {
        if (img->pixDims[i]!=maskImg->pixDims[i]) {
            msg_error("Mask and template do not have same voxel dimensions.");
            return false;
        }

        if (img->imgDims[i]!=maskImg->imgDims[i]) {
            msg_error("Mask and template do not have same image dimensions.");
            return false;
        }

        for (int j=0; j<4; j++)
            if (img->xyz2ijk[i][j]!=maskImg->xyz2ijk[i][j]) {
                msg_error("Mask and template images do not have the same coordinate space.");
                return false;
            }
    }

    maskImg->read();
    
    mask = new bool**[maskImg->imgDims[0]];
    for (int i = 0; i < maskImg->imgDims[0]; i++) {
        mask[i] = new bool*[maskImg->imgDims[1]];
        for (int j = 0; j < maskImg->imgDims[1]; j++) {
            mask[i][j] = new bool[maskImg->imgDims[2]];
            for (int k = 0; k < maskImg->imgDims[2]; k++)
                mask[i][j][k] = ((*maskImg)(i,j,k)>0) ? true : false;
        }
    }

    return true;

}


template<typename T>
bool Tractogram2ImageMapper<T>::setMask(Image<int>* maskImg, int selectedLabel) {
        
    maskImg->readHeader();

    for (int i=0; i<3; i++) {
        if (img->pixDims[i]!=maskImg->pixDims[i]) {
            msg_error("Mask and template do not have same voxel dimensions.");
            return false;
        }

        if (img->imgDims[i]!=maskImg->imgDims[i]) {
            msg_error("Mask and template do not have same image dimensions.");
            return false;
        }

        for (int j=0; j<4; j++)
            if (img->xyz2ijk[i][j]!=maskImg->xyz2ijk[i][j]) {
                msg_error("Mask and template images do not have the same coordinate space.");
                return false;
            }
    }

    maskImg->read();

    mask = new bool**[maskImg->imgDims[0]];
    for (int i = 0; i < maskImg->imgDims[0]; i++) {
        mask[i] = new bool*[maskImg->imgDims[1]];
        for (int j = 0; j < maskImg->imgDims[1]; j++) {
            mask[i][j] = new bool[maskImg->imgDims[2]];
            for (int k = 0; k < maskImg->imgDims[2]; k++)
                mask[i][j][k] = ((*maskImg)(i,j,k)==selectedLabel) ? true : false;
        }
    }

    return true;

}


template<typename T>
void Tractogram2ImageMapper<T>::run(
        std::function<void(Tractogram2ImageMapper<T>* tim, uint16_t _threadNo, int* _gridPos, Segment& _seg)> processor_f,
        std::function<void(Tractogram2ImageMapper<T>* tim)> outputCompiler_f
        ) {

    // Process the tractogram and fill the grid
    if (!QUITE)
    MT::MTRUN(tractogram[0].numberOfStreamlines, MT::maxNumberOfThreads, "Mapping tractogram", [&](MTTASK task)->void {
        processStreamline(task.no,task.threadId, processor_f);
        });
    else
    MT::MTRUN(tractogram[0].numberOfStreamlines, MT::maxNumberOfThreads, [&](MTTASK task)->void {
        processStreamline(task.no,task.threadId, processor_f);
        });

    // Compile output
    outputCompiler_f(this);
    
}

template<typename T1>
bool Tractogram2ImageMapper<T1>::processStreamline(int streamlineId, uint16_t threadNo, std::function<void(Tractogram2ImageMapper<T1>* tim, uint16_t threadNo, int* gridPos, Segment& seg)> f ) {
    
    // If streamline is empty, exit.
    int len = tractogram[threadNo].len[streamlineId];

    if (len==0) return true;
    
    float p0[3], p1[3], dir[3], length, lengthR, lengthScale, t;

    int32_t A[3], B[3];

    // Make streamline kernel, which is used if anisotropic smoothing is applied
    std::vector<float**> kernel;
    kernel.resize(std::get<1>(smoothing)+1);
    kernel[0] = tractogram[threadNo].readStreamline(streamlineId);
    
    if (std::get<0>(smoothing) > 0 ) {

        // Allocate memory for the kernel
        for (int n=0; n<std::get<1>(smoothing); n++) {
            kernel[n+1] = new float*[len];
            for (int l=0; l<len; l++)
                kernel[n+1][l] = new float[3];
        }

        // Prepare points to track
        std::vector<float> scale_N;
        std::vector<float> scale_B;
        for (int n=0; n<std::get<1>(smoothing); n++) {
            scale_N.push_back(MT::rndm[threadNo]->normal_m0_s1()*std::get<0>(smoothing));
            scale_B.push_back(MT::rndm[threadNo]->normal_m0_s1()*std::get<0>(smoothing));
        }

        // To make the kernel, we will compute the parallel transport frame
        float T[3], N[3], B[3], curr_T[3];
        float proj, R_axis[3], R_angle;
        float R[4][4];

        // Get a random initial PTF and handle the initial point
        vec3sub(T,kernel[0][1],kernel[0][0]);
        normalize(T);
        MT::rndm[threadNo]->getAUnitRandomPerpVector(N,T);
        normalize(N);
        cross(B,T,N);

        for (int n=0; n<std::get<1>(smoothing); n++)
            for (int i=0; i<3; i++)
                kernel[n+1][0][i] = kernel[0][0][i] + N[i]*scale_N[n] + + B[i]*scale_B[n];

        for (auto l=1; l<(len-1); l++) {

            vec3sub(curr_T,kernel[0][l+1],kernel[0][l]);
            normalize(curr_T);

            proj       = dot(T,curr_T);
            if (proj >  1) proj =  1;
            if (proj < -1) proj = -1;
            R_angle = std::acos(proj);

            curr_T[0] += EPS6;  // for numerical stability reasons
            curr_T[1] += EPS6;
            curr_T[2] += EPS6;
            
            cross(R_axis,T,curr_T);
            normalize(R_axis);

            axisangle2Rotation(R_axis, R_angle, R);
            rotate(curr_T,T,R); // curr_T is used as a temp var
            rotate(R_axis,N,R); // R_axis is used as a temp var

            T[0] = curr_T[0];
            T[1] = curr_T[1];
            T[2] = curr_T[2];
            normalize(T);

            N[0] = R_axis[0];
            N[1] = R_axis[1];
            N[2] = R_axis[2];
            normalize(N);

            cross(B,T,N);

            for (int n=0; n<std::get<1>(smoothing); n++)
                for (int i=0; i<3; i++)
                    kernel[n+1][l][i] = kernel[0][l][i] + N[i]*scale_N[n] + + B[i]*scale_B[n];

        }

        // Repeat the last N, B for the last segment
        for (int n=0; n<std::get<1>(smoothing); n++)
            for (int i=0; i<3; i++)
                kernel[n+1][len-1][i] = kernel[0][len-1][i] + N[i]*scale_N[n] + + B[i]*scale_B[n];

        
    }

    // Iterate over the kernel
    for (float** streamline : kernel) {
    
        Segment seg;
        seg.streamlineNo = streamlineId;
        
        // End of segment in real and its corner in image space
        img->to_ijk(streamline[0],p0);
        A[0]  = std::round(p0[0]);
        A[1]  = std::round(p0[1]);
        A[2]  = std::round(p0[2]);
        
        // If streamline has 1 point and no segment
        if (len==1) {
            if ( img->isInside(A) && ((mask==NULL) || mask[A[0]][A[1]][A[2]]) ) {
                
                for (int m=0;m<3;m++)
                    seg.p[m] = streamline[0][m];
                seg.dir[0] = 1;
                seg.dir[1] = 0;
                seg.dir[2] = 0;
                seg.length = 0;
                
                f(this, threadNo, A, seg);
                
            }
            delete[] streamline[0];
            delete[] streamline;
            continue;
        }
        
        auto sub2ind = [&] () {return A[0] + A[1]*img->imgDims[1] + A[2]*img->imgDims[1]*img->imgDims[2];};
        
        std::set<size_t> voxelInds;
        std::pair<std::set<size_t>::iterator,bool> voxelChecker;
        
        // If streamline has many points and segments
        for (int i=0; i<len-1; i++) {

            // End of segment in real and its corner in image space
            img->to_ijk(streamline[i+1],p1);
            for (int m=0;m<3;m++) {
                seg.p[m]    = streamline[i][m];
                seg.dir[m]  = streamline[i+1][m] - streamline[i][m];
                B[m]        = std::round(p1[m]);
            }    
            
            // Find segment length and direction in real space
            lengthR = norm(seg.dir);
            for (int m=0;m<3;m++)
                seg.dir[m] /= lengthR;
            
            
            // Segment does not leave the voxel, add this segment and continue with the next one
            if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {
                
                if ( img->isInside(A) && ((mask==NULL) || mask[A[0]][A[1]][A[2]]) ) {
                    
                    seg.length = lengthR;
                    
                    if (mapOnce) {
                        voxelChecker = voxelInds.insert(sub2ind());
                        if (voxelChecker.second) {
                            f(this, threadNo, A, seg);
                        }
                    } else {
                        f(this, threadNo, A, seg);
                    }
                        
                    
                }
                
                for (int m=0;m<3;m++)
                    p0[m]      = p1[m];
                
                continue;
                
            }
            
            // Find length and direction of segment in grid space
            for (int m=0;m<3;m++)
                dir[m]     = p1[m] - p0[m];
            length  = norm(dir);
            for (int m=0;m<3;m++)
                dir[m]     /= length;
            
            
            // Grid lengthScale
            lengthScale = std::sqrt(dir[0]*img->pixDims[0]*dir[0]*img->pixDims[0]+dir[1]*img->pixDims[1]*dir[1]*img->pixDims[1]+dir[2]*img->pixDims[2]*dir[2]*img->pixDims[2]);
            
            t = 0;
            
            while (1) {
            
                if (t>0) {
                    
                    t += EPS4;
                    
                    if ( img->isInside(A) && ((mask==NULL) || mask[A[0]][A[1]][A[2]]) ) {
                        
                        seg.length = t*lengthScale;
                        
                        if (lengthR<seg.length) {
                            seg.length = lengthR;
                            lengthR    = -1;
                        } else {
                            lengthR   -= seg.length;
                        }
                        
                        if (mapOnce) {
                            voxelChecker = voxelInds.insert(sub2ind());
                            if (voxelChecker.second) {
                                f(this, threadNo, A, seg);
                            }
                        } else {
                            f(this, threadNo, A, seg);
                        }
                        
                        for (int m=0;m<3;m++)
                            seg.p[m] += seg.length*seg.dir[m];
                    }

                    
                    length -= t;
                    if ((length<0) || (lengthR<0))
                        break;
                
                    for (int m=0;m<3;m++) {
                        p0[m] += t*dir[m];
                        A[m]   = std::round(p0[m]);
                    }
                    
                    // Remaining part does not leave the voxel, add this part and continue with the next segment
                    if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {
                        if ( img->isInside(A) && ((mask==NULL) || mask[A[0]][A[1]][A[2]]) ) {
                            seg.length = lengthR;
                            if (mapOnce) {
                                voxelChecker = voxelInds.insert(sub2ind());
                                if (voxelChecker.second) {
                                    f(this, threadNo, A, seg);
                                }
                            } else {
                                f(this, threadNo, A, seg);
                            }
                        }
                        break;
                    }
                    
                }
                
                
                if (rayTrace_box(A,p0,dir,t)==false)
                    continue;

                if (length<EPS4)
                    break;
                
                length -= EPS4;
                
                for (int m=0;m<3;m++) {
                    p0[m] += EPS4*dir[m];
                    A[m]   = std::round(p0[m]);
                }
                
                t = 0;
            
            }

            for (int m=0;m<3;m++) {
                p0[m] = p1[m];
                A[m]  = B[m];
            }
            

        }
        
        
        for (int l=0; l<len; l++)
            delete[] streamline[l];
        delete[] streamline;

    }
    
    return true;
}
