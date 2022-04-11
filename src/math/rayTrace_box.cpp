#include "core.h"

bool rayTrace_box(int* A, float* p0, float* dir, float &t) {

    float T[3],u,v;
    
    for (int m=0;m<3;m++)
        T[m] = p0[m] - A[m] + 0.5f;
    
    if (dir[0]<0.0f) {
        // Face 1
        t = -T[0]/dir[0];
        if (t>0.0f) {
            u = T[1]+t*dir[1];
            v = T[2]+t*dir[2];
            if ( ((u<0.0f)||(u>1.0f)||(v<0.0f)||(v>1.0f)) == false)
                return false;
        }
    } else {
        // Face 2
        t = -(T[0]-1)/dir[0];
        if (t>0.0f) {
            u = T[2]+t*dir[2];
            v = T[1]+t*dir[1];
            if ( ((u<0.0f)||(u>1.0f)||(v<0.0f)||(v>1.0f)) == false)
                return false;
        }
    }
    
    
    if (dir[2]<0.0f) {
        // Face 3
        t = -T[2]/dir[2];
        if (t>0.0f) {
            u = T[0]+t*dir[0];
            v = T[1]+t*dir[1];
            if ( ((u<0.0f)||(u>1.0f)||(v<0.0f)||(v>1.0f)) == false)
                return false;
        }
    } else {
        // Face 4
        t = -(T[2]-1)/dir[2];
        if (t>0.0f) {
            u = T[1]+t*dir[1];
            v = T[0]+t*dir[0];
            if ( ((u<0.0f)||(u>1.0f)||(v<0.0f)||(v>1.0f)) == false)
                return false;
        }
    }
    
    
    if (dir[1]<0.0f) {
        // Face 5
        t = -T[1]/dir[1];
        if (t>0.0f) {
            u = T[2]+t*dir[2];
            v = T[0]+t*dir[0];
            if ( ((u<0.0f)||(u>1.0f)||(v<0.0f)||(v>1.0f)) == false)
                return false;
        }
    } else {
        // Face 6
        t = -(T[1]-1)/dir[1];
        if (t>0.0f) {
            u = T[0]+t*dir[0];
            v = T[2]+t*dir[2];
            if ( ((u<0.0f)||(u>1.0f)||(v<0.0f)||(v>1.0f)) == false)
                return false;
        }
    }
    
    t = 0.0f;
    return true;
    
}