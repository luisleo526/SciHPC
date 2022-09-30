//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_GLOBAL_H
#define SCIHPC_GLOBAL_H

typedef double DataType;

struct indices {
    int i, j, k;
};

struct axis {
    DataType start, end;
    int n;
};

DataType ***init_array(int Nx, int Ny, int Nz);

void delete4d(DataType ****arr, int Nx, int Ny, int Nz);

void delete3d(DataType ***arr, int Nx, int Ny);

void delete2d(DataType **arr, int Nx);

#endif //SCIHPC_GLOBAL_H
