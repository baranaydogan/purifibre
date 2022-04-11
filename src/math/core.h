#ifndef SRC_MATH_CORE_H_
#define SRC_MATH_CORE_H_


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <cfloat>
#include <functional>


#define SQRT2	   1.41421356237309504880168872420969807856967187537694807
#define PI 		   3.14159265358979323846264338327950288419716939937510582
#define TWOPI	   6.28318530717958647692528676655900576839433879875021164
#define PIOVERTWO  1.57079632679489661923132169163975144209858469968755291
#define N180OVERPI 57.2957795130823208767981548141051703324054724665643215 // 180/PI
#define PIOVERN180 0.01745329251994329576923690768488612713442871888541725 // PI/180

#define EPS8       0.00000001
#define EPS7       0.0000001
#define EPS6       0.000001
#define EPS5       0.00001
#define EPS4       0.0001
#define EPS3       0.001
#define EPS2       0.01


typedef struct {
  float x;
  float y;
  float z;
} Point;

template<class T>
inline float norm(const T v) 
{ 
    return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

template<class T>
inline void normalize(T v) 
{
    float scale = 1.0/std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] *= scale;
    v[1] *= scale;
    v[2] *= scale;
}



template<class T1,class T2>
inline float dot(const T1 v1,const T2 v2) 
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}


template<class TOUT,class T1,class T2>
inline void cross(TOUT out, const T1 v1,const T2 v2)
{
    out[0] = v1[1]*v2[2] - v1[2]*v2[1];
    out[1] = v1[2]*v2[0] - v1[0]*v2[2];
    out[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

template<class T1,class T2>
inline float dist(const T1 p1,const T2 p2)
{
    return std::sqrt( (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]) );
}

inline float dist(const Point* p1,const Point* p2)
{
    return std::sqrt( (p1->x-p2->x)*(p1->x-p2->x) + (p1->y-p2->y)*(p1->y-p2->y) + (p1->z-p2->z)*(p1->z-p2->z) );
}

inline float dist(const Point p1,const Point p2)
{
    return std::sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z) );
}



template<class TOUT,class T1,class T2>
inline void vec3sub(TOUT out, const T1 v1,const T2 v2)
{
    out[0] = v1[0] - v2[0];
    out[1] = v1[1] - v2[1];
    out[2] = v1[2] - v2[2];
}

template<class TOUT,class T1,class T2>
inline void vec3add(TOUT out, const T1 v1,const T2 v2)
{
    out[0] = v1[0] + v2[0];
    out[1] = v1[1] + v2[1];
    out[2] = v1[2] + v2[2];
}

template<class TOUT,class T1,class T2, class T3>
inline void vec3add(TOUT out, const T1 v1,const T2 v2, const T3 s)
{
    out[0] = v1[0] + s*v2[0];
    out[1] = v1[1] + s*v2[1];
    out[2] = v1[2] + s*v2[2];
}

template<class T1,class T2>
inline void vec3scale(T1 v, const T2 s)
{
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
}

template<class TOUT>
inline void point2vec3(TOUT out, Point p, Point ref)
{
    out[0] = p.x - ref.x;
    out[1] = p.y - ref.y;
    out[2] = p.z - ref.z;
}


template<class T>
inline void applyTransform(Point p, T M)
{
    float x = p.x*M[0][0] + p.y*M[0][1] + p.z*M[0][2] + M[0][3];
    float y = p.x*M[1][0] + p.y*M[1][1] + p.z*M[1][2] + M[1][3];
    float z = p.x*M[2][0] + p.y*M[2][1] + p.z*M[2][2] + M[2][3];
    p.x = x;
    p.y = y;
    p.z = z;
}

template<class T1,class T2>
inline void applyTransform(T1 v, T2 M)
{
    float x = v[0]*M[0][0] + v[1]*M[0][1] + v[2]*M[0][2] + M[0][3];
    float y = v[0]*M[1][0] + v[1]*M[1][1] + v[2]*M[1][2] + M[1][3];
    float z = v[0]*M[2][0] + v[1]*M[2][1] + v[2]*M[2][2] + M[2][3];
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

template<class T1,class T2>
inline void applyTransform(std::vector<T1>& v_arr, T2 M)
{

    for (auto v : v_arr) {
        float x = v[0]*M[0][0] + v[1]*M[0][1] + v[2]*M[0][2] + M[0][3];
        float y = v[0]*M[1][0] + v[1]*M[1][1] + v[2]*M[1][2] + M[1][3];
        float z = v[0]*M[2][0] + v[1]*M[2][1] + v[2]*M[2][2] + M[2][3];
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }
    
}


template<class T>
T rad2deg(T rad) {
    return rad*N180OVERPI;
}

template<class T>
T deg2rad(T deg) {
    return deg*PIOVERN180;
}

template<class T>
void disp(T v)
{
    std::cout << "[ " << v[0] << " , " << v[1] << " , " << v[2] << " ]" << std::endl << std::flush;
}




// All is done in image space
// A: voxel indices, i.e. center of voxel in image space
// p0: origin of ray
// dir: direction of ray
// t: if ray intersects the box, t is the distance from the ray origin
// returns false if there is intersection. 
bool rayTrace_box(int* A, float* p0, float*  dir, float  &t);

#endif
