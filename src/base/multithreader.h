#ifndef SRC_MULTITHREADER_H_
#define SRC_MULTITHREADER_H_

#ifdef BUILD_FOR_WINDOWS
#include <windows.h>
#include <io.h>
#undef max
#else
#include <unistd.h>
#endif


#include <iostream>
#include <string>
#include <iomanip>
#include <thread>
#include <vector>
#include <memory>
#include <functional> 
#include <mutex>
#include <condition_variable>
#include "../math/randomThings.h"

struct MTTASK {
    size_t   no;
    uint16_t threadId;
};

namespace MT {
    
extern std::condition_variable                   exit_cv;
extern std::mutex                                exit_mx;
extern std::mutex                                proc_mx;
extern float                                     progress;
extern int                                       maxNumberOfThreads;
extern uint16_t                                  finishedThreadId;
extern std::vector<std::unique_ptr<RandomDoer>>  rndm;

extern void MTINIT();

extern void MTRUN(size_t range, int numberOfThreads, std::string message, std::function<void(MTTASK mttask)> f);
extern void MTRUN(size_t range, int numberOfThreads, std::function<void(MTTASK mttask)> f);

extern void MTRUN(size_t range, std::string message, std::function<void(MTTASK mttask)> f);
extern void MTRUN(size_t range, std::function<void(MTTASK mttask)> f);


}

#endif
