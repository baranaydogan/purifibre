#include "config.h"
#include "multithreader.h"

using namespace MT;

std::string SGNTR  = "purifibre v0.1";
std::string FOOTER = "If you use purifibre in your work, don't forget to cite:\nAydogan D. B. \"Fiber coupling (FICO) measure using anisotropic smoothing of track orientation density images for tractogram filtering\", ISMRM 2022\n\n©Copyright 2022, Dogu Baran Aydogan, baran.aydogan@uef.fi";

bool        QUITE  = false;

void INIT() {
    MT::MTINIT();
}
