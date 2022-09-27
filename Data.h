//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_DATA_H
#define SCIHPC_DATA_H

#include "global.h"

class Data {

private:
    int ghc = 3;
    int _Nx, _Ny, _Nz;
public:

    Data(const int nx);

    Data(const int nx, const int ny);

    Data(const int nx, const int ny, const int nz);

    int Nx();

    int Ny();

    int Nz();

    DataType ***data;
    int ndim;

    MapIndex index_mapping(int i, int j, int k);

};


#endif //SCIHPC_DATA_H
