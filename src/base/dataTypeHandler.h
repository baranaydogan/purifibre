#ifndef SRC_DATATYPEHANDLER_H_
#define SRC_DATATYPEHANDLER_H_

#include <iostream>
#include <typeinfo>
#include <unordered_map>
#include <string>
#include <memory>
#include <complex>

typedef enum {
    UNKNOWNDATATYPE,
    BOOL,
    INT8,
    UINT8,
    INT16,
    UINT16,
    INT32,
    UINT32,
    INT64,
    UINT64,
    FLOAT32,
    FLOAT64,
    FLOAT128,
    COMPLEX64,
    COMPLEX128,
    COMPLEX256
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
