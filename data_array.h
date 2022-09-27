//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_DATA_ARRAY_H
#define SCIHPC_DATA_ARRAY_H

#include "global.h"

class data_array {
public:
    int ghc = 3;
    int Nx, Ny, Nz;
    int nx, ny, nz;
    int ndim;
    DataType ***data;
    indices dis;

    data_array(const int _nx);

    data_array(const int _nx, const int _ny);

    data_array(const int _nx, const int _ny, const int _nz);


    indices index_mapping(int i, int j, int k);

};


#endif //SCIHPC_DATA_ARRAY_H
