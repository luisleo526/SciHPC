//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_SCALAR_DATA_H
#define SCIHPC_SCALAR_DATA_H

#include "global.h"

class scalar_data {
public:
    int ghc = 3;
    int Nx, Ny, Nz;
    int nx, ny, nz;
    int ndim;

    DataType ***data, ***flux, ***old;
    DataType ***fx, ***fy, ***fz;
    DataType ***fxx, ***fyy, ***fzz;
    DataType ***fxy, ***fxyy;
    DataType ***fyz, ***fyzz;
    DataType ***fzx, ***fzxx;

    explicit scalar_data(int _nx);
    scalar_data(int _nx, int _ny);
    scalar_data(int _nx, int _ny, int _nz);
    indices index_mapping(int i, int j, int k);
    void store();

};


#endif //SCIHPC_SCALAR_DATA_H
