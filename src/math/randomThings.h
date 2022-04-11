#ifndef SRC_MATH_RANDOMTHINGS_H_
#define SRC_MATH_RANDOMTHINGS_H_

#include <cfloat>
#include <random>

#include "./core.h"
#include "./rotation.h"

class RandomDoer {

public:
	RandomDoer();
	~RandomDoer();

	float 		 uniform_01();		// random number between  0   and 1    uniformly distributed
	float 		 uniform_m05_p05(); // random number between -0.5 and 0.5  uniformly distributed
	float 		 uniform_m1_p1();   // random number between -1   and 1    uniformly distributed
	float 		 normal_m0_s1();    // normal distribution with mean=0 and standard deviation=1

	void         init_uniform_int(int limit);
	int          uniform_int();

	void 		 getAUnitRandomVector(float* out);
	void 		 getAUnitRandomPerpVector(float* out, float* inp);
	void 		 getAUnitRandomQuaternion(float* q);
	void         getARandomRotationMatrix(float R[][4]);
    void         getARandomPointWithinDisk(float* x, float *y, float r);
    void         getARandomPointWithinSphere(float* x, float *y, float *z, float r);
	void 		 getARandomPointWithinTriangle(float* out, float* a, float* b, float* c);

	void 		 getARandomPointWithinBoundingBox(float* p, float* bb); // bb must have values in order, e.g., [x_min x_max y_min y_max z_min z_max] 
	void 		 getNRandomPointsWithinBoundingBox(std::vector<float*>& p, float* bb, int N); // bb must have values in order, e.g., [x_min x_max y_min y_max z_min z_max] 

	std::mt19937 getGen() { return gen; }

private:
	std::mt19937 gen;
	std::uniform_real_distribution<float> *unidis_01;
	std::uniform_real_distribution<float> *unidis_m05_p05;
	std::uniform_real_distribution<float> *unidis_m1_p1;
	std::uniform_int_distribution<int>    *unidis_int;
	std::normal_distribution<float>       *normdis_m0_s1;

};

inline float RandomDoer::uniform_01() {
	return (*unidis_01)(gen);
}

inline float RandomDoer::uniform_m05_p05() {
	return (*unidis_m05_p05)(gen);
}

inline float RandomDoer::uniform_m1_p1() {
	return (*unidis_m1_p1)(gen);
}

inline float RandomDoer::normal_m0_s1() {
	return (*normdis_m0_s1)(gen);
}

inline int RandomDoer::uniform_int() {
	return (*unidis_int)(gen);
}

inline void RandomDoer::getAUnitRandomVector(float* out) {

	do {
		out[0] = (*unidis_m1_p1)(gen);
		out[1] = (*unidis_m1_p1)(gen);
		out[2] = (*unidis_m1_p1)(gen);
	} while ((out[0]==0) && (out[1]==0) && (out[2]==0));
	normalize(out);

}

inline void RandomDoer::getARandomPointWithinDisk(float* x, float *y, float r) {

	do {
		*x = (*unidis_m1_p1)(gen);
		*y = (*unidis_m1_p1)(gen);
	} while ((*x**x+*y**y)>1);
	*x *= r;  
    *y *= r;
}

inline void RandomDoer::getARandomPointWithinSphere(float* x, float *y, float *z, float r) {

	do {
		*x = (*unidis_m1_p1)(gen);
		*y = (*unidis_m1_p1)(gen);
        *z = (*unidis_m1_p1)(gen);
	} while ((*x**x+*y**y+*z**z)>1);
	*x *= r;  
    *y *= r;
    *z *= r;
}

inline void RandomDoer::getARandomPointWithinTriangle(float* out, float* a, float* b, float* c) {

	float r1 = (*unidis_01)(gen);
	float r2 = (*unidis_01)(gen);

	// Generate a uniformly random point in the triangle with vertices (0,0) (0,1) (1,0).
	// If the point (r1, r2) is not inside the triangle, reflect it through the point (0.5, 0.5).
	if (r1 + r2 > 1) {
		r1 = 1 - r1;
		r2 = 1 - r2;
	}

	// Map the point (r1, r2) in the unit-leg triangle to a point in the triangle ABC.
	// Intuition: the unit legs are stretched to match the sides AB and AC, and the origin
	// is moved to A.
	out[0] = a[0] + r1 * (b[0]-a[0]) + r2 * (c[0]-a[0]);
	out[1] = a[1] + r1 * (b[1]-a[1]) + r2 * (c[1]-a[1]);
	out[2] = a[2] + r1 * (b[2]-a[2]) + r2 * (c[2]-a[2]);
}

inline void RandomDoer::getAUnitRandomPerpVector(float* out, float* inp) {
	float tmp[3];
	getAUnitRandomVector(tmp);
	cross(out,inp,tmp);
	normalize(out);
}

inline void RandomDoer::getAUnitRandomQuaternion(float* q) {
	q[0] = (*normdis_m0_s1)(gen);
	q[1] = (*normdis_m0_s1)(gen);
	q[2] = (*normdis_m0_s1)(gen);
	q[3] = (*normdis_m0_s1)(gen);
	float normalizer = 1.0/std::sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
	q[0] *= normalizer;
	q[1] *= normalizer;
	q[2] *= normalizer;
	q[3] *= normalizer;
}

inline void RandomDoer::getARandomRotationMatrix(float R[][4]) {
	float q[4];
	getAUnitRandomQuaternion(&q[0]);
	quaternion2Rotation(&q[0],R);
}

// bb must have values in order, e.g., [x_min x_max y_min y_max z_min z_max] 
inline void RandomDoer::getARandomPointWithinBoundingBox(float* p, float* bb) {
	p[0] = (*unidis_01)(gen)*(bb[1]-bb[0]) + bb[0];
	p[1] = (*unidis_01)(gen)*(bb[3]-bb[2]) + bb[2];
	p[2] = (*unidis_01)(gen)*(bb[5]-bb[4]) + bb[4];
}

// bb must have values in order, e.g., [x_min x_max y_min y_max z_min z_max] 
inline void RandomDoer::getNRandomPointsWithinBoundingBox(std::vector<float*>& p, float* bb, int N) {

	p.reserve(N);

	for (size_t i=0; i<p.size(); i++) {
		float* r = new float[3];
		getARandomPointWithinBoundingBox(r,bb);
		p.push_back(r);
	}

}

#endif
