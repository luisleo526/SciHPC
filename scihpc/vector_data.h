//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_VECTOR_DATA_H
#define SCIHPC_VECTOR_DATA_H

#include "scalar_data.h"

class vector_data {
public:
    vector_data(const int _nx);

    vector_data(const int _nx, const int _ny);

    vector_data(const int _nx, const int _ny, const int _nz);

    scalar_data x, y, z;
};


#endif //SCIHPC_VECTOR_DATA_H