#ifndef SRC_BYTESWAPPER_H_
#define SRC_BYTESWAPPER_H_

template<typename T>
inline void swapByteOrder(T& a) {
    char* byteArray = reinterpret_cast<char*>(&a);
    for(size_t i = 0; i < (sizeof(a)/2); i++)
        std::swap(byteArray[sizeof(a) - 1 - i],byteArray[i]);
}

#endif
