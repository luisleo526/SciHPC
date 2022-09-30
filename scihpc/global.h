//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_GLOBAL_H
#define SCIHPC_GLOBAL_H

typedef long double DataType;

struct indices {
    int i, j, k;
};

struct axis {
    DataType start, end;
    int n;
};

DataType ***init_array(int Nx, int Ny, int Nz);

void delete3d(DataType ***arr, int Nx, int Ny);

void delete2d(DataType **arr, int Nx);

DataType* x_at(int Nx, int j, int k, DataType ***arr);
DataType* y_at(int Ny, int i, int k, DataType ***arr);
DataType* z_at(int Nz, int i, int j, DataType ***arr);

#endif //SCIHPC_GLOBAL_H
