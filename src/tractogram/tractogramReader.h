#ifndef SRC_TRACTOGRAMREADER_H_
#define SRC_TRACTOGRAMREADER_H_

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <float.h>
#include "../base/config.h"
#include "../math/core.h"

struct Segment {
    int     streamlineNo;
    float   p[3];
    float   dir[3];
    float   length;
    void*   data;
};

typedef enum {
    VTK_ASCII,
    VTK_BINARY,
    TCK,
    TRK
} TRACTOGRAMFILEFORMAT;

class TractogramReader {
    
public:
    
    TractogramReader();
    TractogramReader(std::string _fileName);
    ~TractogramReader();
    TractogramReader(const TractogramReader& obj);
    void copyFrom(const TractogramReader& obj);
    void destroyCopy();
    
    bool     initReader(std::string _fileName);
    float**  readStreamline(size_t n);
    void     readPoint(size_t n, uint32_t l, float* p);
    void     setThreadId (uint32_t _threadId) {threadId=_threadId;}
    uint32_t getThreadId () {return threadId;}
    
    std::vector<Point> readVectorizedStreamline(size_t n);

    FILE                   *file;
    std::string             fileName;
    std::string             fileDescription;
    TRACTOGRAMFILEFORMAT    fileFormat;
    
    size_t                  numberOfPoints;
    size_t                  numberOfStreamlines;
    uint32_t*               len;
    
    
private:
    
    long*                   streamlinePos;  // file positions for first points of streamlines
    uint32_t                threadId;

};

#endif

