#ifndef SRC_VECTOROPERATIONS_H_
#define SRC_VECTOROPERATIONS_H_

#include <algorithm>
#include <vector>

template<class T1,class T2>
void removeIdx(std::vector<T1>& inp, const std::vector<T2>& rmIdx) {
    
    if (rmIdx.size()==0) {
        return;
    } else if (inp.size()==rmIdx.size()) {
        std::vector<T1>().swap(inp);
        return;
    } 

    std::vector<T2> idx = rmIdx;

    std::sort(idx.begin(),idx.end());

    auto   beg    = inp.begin();
    size_t offset = 0;

    for (auto it = idx.begin(); it < idx.end(); it++) {
        size_t next = (it + 1 == idx.cend() ? inp.size() : *(it + 1));
        std::move(beg+*it+1, beg+next, beg+*it-offset);
        offset++;
    }
    
    inp.resize(inp.size()-idx.size());

}

#endif