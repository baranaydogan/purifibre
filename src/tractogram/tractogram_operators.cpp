#include "tractogram_operators.h"

std::vector<float> getTractogramBBox(TractogramReader* _tractogram) {

    std::vector<float> bb(6,0);

    std::vector<std::vector<float>> thread_bb; 

    // Initialize tractogram and make copies for multithreader
    TractogramReader* tractogram = new TractogramReader[MT::maxNumberOfThreads]();
    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        tractogram[t].copyFrom(*_tractogram);
        std::vector<float> tmp(6,0);
        thread_bb.push_back(tmp);
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

    for (int t = 0; t < MT::maxNumberOfThreads; t++)
        for (int i=0; i<6; i++)
            thread_bb[t][i] = bb[i];
    

    for (uint32_t i=0; i<tractogram[0].len[0]; i++)
        delete[] firstStreamline[i];
    delete[] firstStreamline;

    // Iterate throught the whole tractogram
    auto findBB = [&](MTTASK task)->void{
        
        float** streamline = tractogram[task.threadId].readStreamline(task.no);

        for (uint32_t i=0; i<tractogram[task.threadId].len[task.no]; i++) {

            if ( streamline[i][0] < thread_bb[task.threadId][0] )   thread_bb[task.threadId][0] = streamline[i][0];
            if ( streamline[i][0] > thread_bb[task.threadId][1] )   thread_bb[task.threadId][1] = streamline[i][0];
            if ( streamline[i][1] < thread_bb[task.threadId][2] )   thread_bb[task.threadId][2] = streamline[i][1];
            if ( streamline[i][1] > thread_bb[task.threadId][3] )   thread_bb[task.threadId][3] = streamline[i][1];
            if ( streamline[i][2] < thread_bb[task.threadId][4] )   thread_bb[task.threadId][4] = streamline[i][2];
            if ( streamline[i][2] > thread_bb[task.threadId][5] )   thread_bb[task.threadId][5] = streamline[i][2];
            
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

        if ( thread_bb[t][0] < bb[0] )   bb[0] = thread_bb[t][0];
        if ( thread_bb[t][1] > bb[1] )   bb[1] = thread_bb[t][1];
        if ( thread_bb[t][2] < bb[2] )   bb[2] = thread_bb[t][2];
        if ( thread_bb[t][3] > bb[3] )   bb[3] = thread_bb[t][3];
        if ( thread_bb[t][4] < bb[4] )   bb[4] = thread_bb[t][4];
        if ( thread_bb[t][5] > bb[5] )   bb[5] = thread_bb[t][5];


    }
    delete[] tractogram;

    return bb;

}
