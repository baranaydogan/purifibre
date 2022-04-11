#include "dataTypeHandler.h"


std::unordered_map<DATATYPE, TypeInfoRef> TYPE = {
    {BOOL,      typeid(bool)                        },
    {INT8,      typeid(int8_t)                      },
    {UINT8,     typeid(uint8_t)                     },
    {INT16,     typeid(int16_t)                     },
    {UINT16,    typeid(uint16_t)                    },
    {INT32,     typeid(int32_t)                     },
    {UINT32,    typeid(uint32_t)                    },
    {INT64,     typeid(int64_t)                     },
    {UINT64,    typeid(uint64_t)                    },
    {FLOAT32,   typeid(float)                       },
    {FLOAT64,   typeid(double)                      },
    {FLOAT128,  typeid(long double)                 },
    {COMPLEX64, typeid(std::complex<float>)         },
    {COMPLEX128,typeid(std::complex<double>)        },
    {COMPLEX256,typeid(std::complex<long double>)   }
};



std::unordered_map<TypeInfoRef, DATATYPE, Hasher, EqualTo> TYPEIDS = {
    {typeid(bool),                          BOOL},
    {typeid(int8_t),                        INT8},
    {typeid(uint8_t),                       UINT8},
    {typeid(int16_t),                       INT16},
    {typeid(uint16_t),                      UINT16},
    {typeid(int32_t),                       INT32},
    {typeid(uint32_t),                      UINT32},
    {typeid(int64_t),                       INT64},
    {typeid(uint64_t),                      UINT64},
    {typeid(float),                         FLOAT32},
    {typeid(double),                        FLOAT64},
    {typeid(long double),                   FLOAT128},
    {typeid(std::complex<float>),           COMPLEX64},
    {typeid(std::complex<double>),          COMPLEX128},
    {typeid(std::complex<long double>),     COMPLEX256}
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
