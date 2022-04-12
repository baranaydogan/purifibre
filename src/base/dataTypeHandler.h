#ifndef SRC_DATATYPEHANDLER_H_
#define SRC_DATATYPEHANDLER_H_

#include <iostream>
#include <typeinfo>
#include <unordered_map>
#include <string>
#include <memory>
#include <complex>

typedef enum {
    UNKNOWN_DT,
    BOOL_DT,
    INT8_DT,
    UINT8_DT,
    INT16_DT,
    UINT16_DT,
    INT32_DT,
    UINT32_DT,
    INT64_DT,
    UINT64_DT,
    FLOAT32_DT,
    FLOAT64_DT,
    FLOAT128_DT,
    COMPLEX64_DT,
    COMPLEX128_DT,
    COMPLEX256_DT
} DATATYPE;


using TypeInfoRef = std::reference_wrapper<const std::type_info>;

struct Hasher {
    std::size_t operator()(TypeInfoRef code) const
    {
        return code.get().hash_code();
    }
};

struct EqualTo {
    bool operator()(TypeInfoRef lhs, TypeInfoRef rhs) const
    {
        return lhs.get() == rhs.get();
    }
};

extern std::unordered_map<DATATYPE, TypeInfoRef>                     TYPE;
extern std::unordered_map<TypeInfoRef, DATATYPE,    Hasher, EqualTo> TYPEIDS;
extern std::unordered_map<TypeInfoRef, std::string, Hasher, EqualTo> TYPENAMES;


#endif
