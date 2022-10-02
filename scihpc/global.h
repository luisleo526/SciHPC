//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_GLOBAL_H
#define SCIHPC_GLOBAL_H

#include <cmath>

typedef double DataType;

const DataType pi = acos(static_cast<DataType>(-1.0));

struct indices {
    int i, j, k;
};

struct axis {
    DataType start, end;
    int n;
};

struct problem_parameters{
    DataType dt;
    DataType density_ratio;
    DataType viscosity_ratio;
    DataType ls_width;
    DataType lsf_mass0;
};

DataType ***init_array(int Nx, int Ny, int Nz);

void delete4d(DataType ****arr, int Nx, int Ny, int Nz);

void delete3d(DataType ***arr, int Nx, int Ny);

void delete2d(DataType **arr, int Nx);

#endif //SCIHPC_GLOBAL_H
