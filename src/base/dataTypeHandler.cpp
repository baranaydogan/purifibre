#include "dataTypeHandler.h"


std::unordered_map<DATATYPE, TypeInfoRef> TYPE = {
    {BOOL_DT,      typeid(bool)                        },
    {INT8_DT,      typeid(int8_t)                      },
    {UINT8_DT,     typeid(uint8_t)                     },
    {INT16_DT,     typeid(int16_t)                     },
    {UINT16_DT,    typeid(uint16_t)                    },
    {INT32_DT,     typeid(int32_t)                     },
    {UINT32_DT,    typeid(uint32_t)                    },
    {INT64_DT,     typeid(int64_t)                     },
    {UINT64_DT,    typeid(uint64_t)                    },
    {FLOAT32_DT,   typeid(float)                       },
    {FLOAT64_DT,   typeid(double)                      },
    {FLOAT128_DT,  typeid(long double)                 },
    {COMPLEX64_DT, typeid(std::complex<float>)         },
    {COMPLEX128_DT,typeid(std::complex<double>)        },
    {COMPLEX256_DT,typeid(std::complex<long double>)   }
};



std::unordered_map<TypeInfoRef, DATATYPE, Hasher, EqualTo> TYPEIDS = {
    {typeid(bool),                          BOOL_DT},
    {typeid(int8_t),                        INT8_DT},
    {typeid(uint8_t),                       UINT8_DT},
    {typeid(int16_t),                       INT16_DT},
    {typeid(uint16_t),                      UINT16_DT},
    {typeid(int32_t),                       INT32_DT},
    {typeid(uint32_t),                      UINT32_DT},
    {typeid(int64_t),                       INT64_DT},
    {typeid(uint64_t),                      UINT64_DT},
    {typeid(float),                         FLOAT32_DT},
    {typeid(double),                        FLOAT64_DT},
    {typeid(long double),                   FLOAT128_DT},
    {typeid(std::complex<float>),           COMPLEX64_DT},
    {typeid(std::complex<double>),          COMPLEX128_DT},
    {typeid(std::complex<long double>),     COMPLEX256_DT}
};

std::unordered_map<TypeInfoRef, std::string, Hasher, EqualTo> TYPENAMES = {
    {typeid(bool),                          "BOOL"},
    {typeid(int8_t),                        "INT8"},
    {typeid(uint8_t),                       "UINT8"},
    {typeid(int16_t),                       "INT16"},
    {typeid(uint16_t),                      "UINT16"},
    {typeid(int32_t),                       "INT32"},
    {typeid(uint32_t),                      "UINT32"},
    {typeid(int64_t),                       "INT64"},
    {typeid(uint64_t),                      "UINT64"},
    {typeid(float),                         "FLOAT32"},
    {typeid(double),                        "FLOAT64"},
    {typeid(long double),                   "FLOAT128"},
    {typeid(std::complex<float>),           "COMPLEX64"},
    {typeid(std::complex<double>),          "COMPLEX128"},
    {typeid(std::complex<long double>),     "COMPLEX256"}
};
