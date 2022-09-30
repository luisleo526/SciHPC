//
// Created by 溫晧良 on 2022/9/28.
//

#include "global.h"

DataType ***init_array(int Nx, int Ny, int Nz) {
    auto ***arr = new DataType **[Nx];
    for (int i = 0; i < Nx; ++i) {
        arr[i] = new DataType *[Ny];
        for (int j = 0; j < Ny; ++j) {
            arr[i][j] = new DataType[Nz];
        }
    }
    return arr;
}

void delete4d(DataType ****arr, int Nx, int Ny, int Nz) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                delete[] arr[i][j][k];
            }
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }
    delete[] arr;
}

void delete3d(DataType ***arr, int Nx, int Ny) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }
    delete[] arr;
}

void delete2d(DataType **arr, int Nx) {
    for (int i = 0; i < Nx; ++i) {
        delete[] arr[i];
    }
    delete[] arr;
}
