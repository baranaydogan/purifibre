#ifndef SRC_IMAGE_INTERPOLATION_H_
#define SRC_IMAGE_INTERPOLATION_H_

template<typename T>
class Image;

typedef enum {
    INSIDE,
    BOUNDARY,
    OUTSIDE
} INTERPAT;

typedef enum {
    NEAREST,
    LINEAR,
    CUBIC
} INTERPMETHOD;

template<typename OUT_T, typename INP_T>
class Interpolator {
    
public:
    
     Interpolator();
    ~Interpolator(){};


    static INTERPAT  init_interp_linear   (Image<INP_T>* img, OUT_T *p, OUT_T* cfs, int64_t* cor_ijk);
    static OUT_T     interp_linear_3D     (Image<INP_T>* img, OUT_T* p);
    static OUT_T     interp_linear_4D_att (Image<INP_T>* img, OUT_T* p, int64_t t);
    static void      interp_linear_4D_all (Image<INP_T>* img, OUT_T* p, OUT_T* out);
    
    
    // TODO: Nearest neighbor interpolation
    static INTERPAT  init_interp_nearest  (Image<INP_T>*,     OUT_T*, OUT_T*, int64_t*) {return OUTSIDE;}
    static OUT_T     interp_nearest_3D    (Image<INP_T>* img, OUT_T*)                   {return img->outsideVal;}
    static OUT_T     interp_nearest_4D_att(Image<INP_T>* img, OUT_T* , int64_t )        {return img->outsideVal;}
    static void      interp_nearest_4D_all(Image<INP_T>*,     OUT_T*, OUT_T*)           {}
    
    // TODO: Cubic interpolation
    static INTERPAT  init_interp_cubic    (Image<INP_T>*,     OUT_T*, OUT_T*, int64_t*) {return OUTSIDE;}
    static OUT_T     interp_cubic_3D      (Image<INP_T>* img, OUT_T*)                   {return img->outsideVal;}
    static OUT_T     interp_cubic_4D_att  (Image<INP_T>* img, OUT_T* , int64_t )        {return img->outsideVal;}
    static void      interp_cubic_4D_all  (Image<INP_T>*,     OUT_T*, OUT_T*)           {}
    
};


// SPEED IS REDUCED IF FUNCTIONS ARE DEFINED IN A SEPARATE CPP FILE
template<typename OUT_T, typename INP_T>
INTERPAT Interpolator<OUT_T,INP_T>::init_interp_linear(Image<INP_T>* img, OUT_T *p, OUT_T* cfs, int64_t* cor_ijk) 
{
    
    OUT_T ijk[3];
    img->to_ijk(p,ijk);
    
    cor_ijk[0]  = int64_t(ijk[0] + 1.0f);
    cor_ijk[1]  = int64_t(ijk[1] + 1.0f);
    cor_ijk[2]  = int64_t(ijk[2] + 1.0f);
    
    cfs[0]      = cor_ijk[0] - ijk[0];
    cfs[1]      = cor_ijk[1] - ijk[1];
    cfs[2]      = cor_ijk[2] - ijk[2];

    if ( cor_ijk[0]>0 && cor_ijk[1]>0 && cor_ijk[2]>0 && cor_ijk[0]<img->imgDims[0] && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]<img->imgDims[2] )
        return INSIDE;
    
    if ( cor_ijk[0]<0 || cor_ijk[1]<0 || cor_ijk[2]<0 || cor_ijk[0]>img->imgDims[0] || cor_ijk[1]>img->imgDims[1] || cor_ijk[2]>img->imgDims[2] )
        return OUTSIDE;
    
    if (cfs[0]>1) cfs[0] -= 1.0f;
    if (cfs[1]>1) cfs[1] -= 1.0f;
    if (cfs[2]>1) cfs[2] -= 1.0f;
    if (cfs[0]<0) cfs[0] += 1.0f;
    if (cfs[1]<0) cfs[1] += 1.0f;
    if (cfs[2]<0) cfs[2] += 1.0f;

    return BOUNDARY;
    
}


template<typename OUT_T, typename INP_T>
OUT_T Interpolator<OUT_T,INP_T>::interp_linear_3D(Image<INP_T>* img, OUT_T* p)
{
    
    OUT_T     cfs[3];
    int64_t   cor_ijk[3];
    
    switch (init_interp_linear(img, p, cfs, cor_ijk)) {
    
        case INSIDE:
            return          cfs[0] *   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1,cor_ijk[2]-1)]) +
                         (1-cfs[0])*   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1,cor_ijk[2]-1)]) +
                            cfs[0] *(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],  cor_ijk[2]-1)]) +
                         (1-cfs[0])*(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],  cor_ijk[2]-1)]) +
                            cfs[0] *   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1,cor_ijk[2]  )]) +
                         (1-cfs[0])*   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1,cor_ijk[2]  )]) +
                            cfs[0] *(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],  cor_ijk[2]  )]) +
                         (1-cfs[0])*(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],  cor_ijk[2]  )]);
            
        case BOUNDARY: {
            
            OUT_T out = 0;
            
            if (cor_ijk[0]>0 && cor_ijk[1]>0 && cor_ijk[2]>0) 
                out += cfs[0] *   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1,cor_ijk[2]-1)]);
            
            if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]>0 && cor_ijk[2]>0) 
                out += (1-cfs[0])*   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1,cor_ijk[2]-1)]);
            
            if (cor_ijk[0]>0 && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]>0)
                out += cfs[0] *(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],  cor_ijk[2]-1)]);
            
            if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]>0)
                out += (1-cfs[0])*(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],  cor_ijk[2]-1)]);
            
            if (cor_ijk[0]>0 && cor_ijk[1]>0 && cor_ijk[2]<img->imgDims[2])
                out += cfs[0] *   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1,cor_ijk[2]  )]);
            
            if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]>0 && cor_ijk[2]<img->imgDims[2])
                out += (1-cfs[0])*   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1,cor_ijk[2]  )]);
            
            if (cor_ijk[0]>0 && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]<img->imgDims[2])
                out += cfs[0] *(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],  cor_ijk[2]  )]);
            
            if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]<img->imgDims[2])
                out += (1-cfs[0])*(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],  cor_ijk[2]  )]);

            return out;
            
            break;
            
            }
        
        default: // OUTSIDE
            return img->outsideVal;
            break;
            
        
    }
    
    return img->outsideVal;

}



template<typename OUT_T, typename INP_T>
OUT_T Interpolator<OUT_T,INP_T>::interp_linear_4D_att(Image<INP_T>* img, OUT_T* p, int64_t t)
{
    OUT_T     cfs[3];
    int64_t   cor_ijk[3];
    
    switch (init_interp_linear(img, p, cfs, cor_ijk)) {
    
        case INSIDE: {

            OUT_T tmp =    cfs[0] *   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1, cor_ijk[2]-1, t)]) +
                        (1-cfs[0])*   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1, cor_ijk[2]-1, t)]) +
                           cfs[0] *(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],   cor_ijk[2]-1, t)]) +
                        (1-cfs[0])*(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],   cor_ijk[2]-1, t)]) +
                           cfs[0] *   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1, cor_ijk[2],   t)]) +
                        (1-cfs[0])*   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1, cor_ijk[2],   t)]) +
                           cfs[0] *(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],   cor_ijk[2],   t)]) +
                        (1-cfs[0])*(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],   cor_ijk[2],   t)]);

            return tmp;

            break;
        }
            
        case BOUNDARY: {
            
            OUT_T tmp = 0;
            
            if (cor_ijk[0]>0 && cor_ijk[1]>0 && cor_ijk[2]>0)
                tmp += cfs[0] *   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1, cor_ijk[2]-1, t)]);
            
            if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]>0 && cor_ijk[2]>0)
                tmp += (1-cfs[0])*   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1, cor_ijk[2]-1, t)]);
            
            if (cor_ijk[0]>0 && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]>0)
                tmp += cfs[0] *(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],   cor_ijk[2]-1, t)]);
            
            if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]>0)
                tmp += (1-cfs[0])*(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],   cor_ijk[2]-1, t)]);
            
            if (cor_ijk[0]>0 && cor_ijk[1]>0 && cor_ijk[2]<img->imgDims[2])
                tmp += cfs[0] *   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1, cor_ijk[2],   t)]);
            
            if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]>0 && cor_ijk[2]<img->imgDims[2])
                tmp += (1-cfs[0])*   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1, cor_ijk[2],   t)]);
            
            if (cor_ijk[0]>0 && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]<img->imgDims[2])
                tmp += cfs[0] *(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],   cor_ijk[2],   t)]);
            
            if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]<img->imgDims[2])
                tmp += (1-cfs[0])*(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],   cor_ijk[2],   t)]);
            
            return tmp;

            break;
            
            }
        
        default: // OUTSIDE
            return img->outsideVal;
            break;
        
    }

    return img->outsideVal;
    
}


template<typename OUT_T, typename INP_T>
void Interpolator<OUT_T,INP_T>::interp_linear_4D_all(Image<INP_T>* img, OUT_T* p, OUT_T* out)
{
    OUT_T     cfs[3];
    int64_t   cor_ijk[3];
    
    switch (init_interp_linear(img, p, cfs, cor_ijk)) {
    
        case INSIDE:
            for (int64_t c=0; c<img->valCnt; c++) {
                out[c] =   cfs[0] *   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1, cor_ijk[2]-1, c)]) +
                        (1-cfs[0])*   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1, cor_ijk[2]-1, c)]) +
                           cfs[0] *(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],   cor_ijk[2]-1, c)]) +
                        (1-cfs[0])*(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],   cor_ijk[2]-1, c)]) +
                           cfs[0] *   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1, cor_ijk[2],   c)]) +
                        (1-cfs[0])*   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1, cor_ijk[2],   c)]) +
                           cfs[0] *(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],   cor_ijk[2],   c)]) +
                        (1-cfs[0])*(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],   cor_ijk[2],   c)]);
            }
            break;
            
        case BOUNDARY: {
            
            for (int64_t c=0; c<img->valCnt; c++) {
            
                OUT_T tmp = 0;
                
                if (cor_ijk[0]>0 && cor_ijk[1]>0 && cor_ijk[2]>0) 
                    tmp += cfs[0] *   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1, cor_ijk[2]-1, c)]);
                
                if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]>0 && cor_ijk[2]>0) 
                    tmp += (1-cfs[0])*   cfs[1] *   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1, cor_ijk[2]-1, c)]);
                
                if (cor_ijk[0]>0 && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]>0)
                    tmp += cfs[0] *(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],   cor_ijk[2]-1, c)]);
                
                if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]>0)
                    tmp += (1-cfs[0])*(1-cfs[1])*   cfs[2] *OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],   cor_ijk[2]-1, c)]);
                
                if (cor_ijk[0]>0 && cor_ijk[1]>0 && cor_ijk[2]<img->imgDims[2])
                    tmp += cfs[0] *   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1]-1, cor_ijk[2],   c)]);
                
                if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]>0 && cor_ijk[2]<img->imgDims[2])
                    tmp += (1-cfs[0])*   cfs[1] *(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1]-1, cor_ijk[2],   c)]);
                
                if (cor_ijk[0]>0 && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]<img->imgDims[2])
                    tmp += cfs[0] *(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0]-1, cor_ijk[1],   cor_ijk[2],   c)]);
                
                if (cor_ijk[0]<img->imgDims[0] && cor_ijk[1]<img->imgDims[1] && cor_ijk[2]<img->imgDims[2])
                    tmp += (1-cfs[0])*(1-cfs[1])*(1-cfs[2])*OUT_T(img->data[img->sub2ind(cor_ijk[0],   cor_ijk[1],   cor_ijk[2],   c)]);
                
                out[c] = tmp;
            }
            
            break;
            
            }
        
        default: // OUTSIDE
            memset(out,img->outsideVal,img->valCnt*sizeof(OUT_T));
            break;
        
    }

    return;
    
}


#endif
