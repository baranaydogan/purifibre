#include "tractogram_operators.h"

std::vector<float> getTractogramBBox(TractogramReader* _tractogram) {

    std::vector<float> bb(6,0); 

    // Initialize tractogram and make copies for multithreader
    TractogramReader* tractogram = new TractogramReader[MT::maxNumberOfThreads]();
    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        tractogram[t].copyFrom(*_tractogram);
    }
    
    int N = tractogram[0].numberOfStreamlines;
    
    if (N<1) {
        return bb;
    }

    // Initialize bb using the first streamline
    float** firstStreamline = tractogram[0].readStreamline(0);
    bb[0] = firstStreamline[0][0];
    bb[1] = firstStreamline[0][0];
    bb[2] = firstStreamline[0][1];
    bb[3] = firstStreamline[0][1];
    bb[4] = firstStreamline[0][2];
    bb[5] = firstStreamline[0][2];

    for (uint32_t i=0; i<tractogram[0].len[0]; i++)
        delete[] firstStreamline[i];
    delete[] firstStreamline;

    // Iterate throught the whole tractogram
    auto findBB = [&](MTTASK task)->void{
        
        float** streamline = tractogram[task.threadId].readStreamline(task.no);

        for (uint32_t i=0; i<tractogram[task.threadId].len[task.no]; i++) {

            if ( streamline[i][0] < bb[0] ) {
                MT::proc_mx.lock();
                bb[0] = streamline[i][0];
                MT::proc_mx.unlock();
            }
            if ( streamline[i][0] > bb[1] ) {
                MT::proc_mx.lock();
                bb[1] = streamline[i][0];
                MT::proc_mx.unlock();
            }
            if ( streamline[i][1] < bb[2] ) {
                MT::proc_mx.lock();
                bb[2] = streamline[i][1];
                MT::proc_mx.unlock();
            }
            if ( streamline[i][1] > bb[3] ) {
                MT::proc_mx.lock();
                bb[3] = streamline[i][1];
                MT::proc_mx.unlock();
            }
            if ( streamline[i][2] < bb[4] ) {
                MT::proc_mx.lock();
                bb[4] = streamline[i][2];
                MT::proc_mx.unlock();
            }
            if ( streamline[i][2] > bb[5] ) {
                MT::proc_mx.lock();
                bb[5] = streamline[i][2];
                MT::proc_mx.unlock();
            }

        }

        for (uint32_t i=0; i<tractogram[task.threadId].len[task.no]; i++)
            delete[] streamline[i];
        delete[] streamline;
        
    };
    if (!QUITE)
        MT::MTRUN(N, MT::maxNumberOfThreads, "Finding tractogram bounding box", findBB);
    else
        MT::MTRUN(N, MT::maxNumberOfThreads, findBB);

    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        tractogram[t].destroyCopy();
    }
    delete[] tractogram;

    return bb;

}
