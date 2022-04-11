#include "randomThings.h"

RandomDoer::RandomDoer() {
	

#ifdef BUILD_FOR_WINDOWS
	std::random_device rd;
	unsigned long randSeed = rd();
#else
	unsigned int lo, hi;
	__asm__ __volatile__("rdtsc" : "=a" (lo), "=d" (hi));
	unsigned long randSeed = ((unsigned long long)hi << 32) | lo;
#endif

	gen.seed(randSeed);
	unidis_01  				= new std::uniform_real_distribution<float>(   0, std::nextafter(1,   FLT_MAX));
	unidis_m05_p05 			= new std::uniform_real_distribution<float>(-0.5, std::nextafter(0.5, FLT_MAX));
	unidis_m1_p1 			= new std::uniform_real_distribution<float>(  -1, std::nextafter(1,   FLT_MAX));
	unidis_int 				= NULL;
	normdis_m0_s1 			= new std::normal_distribution<float>(0.0f,1.0f);

}

RandomDoer::~RandomDoer() {

	delete unidis_01;
	delete unidis_m05_p05;
	delete unidis_m1_p1;
	if (unidis_int!=NULL)
		delete unidis_int;
	delete normdis_m0_s1;
}

void RandomDoer::init_uniform_int(int limit) {
	unidis_int 		= new std::uniform_int_distribution<int>(0,limit);
}
