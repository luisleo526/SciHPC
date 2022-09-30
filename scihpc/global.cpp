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

DataType *x_at(int Nx, int j, int k, DataType ***arr) {
    auto *vector = new DataType[Nx];
    for (int i = 0; i < Nx; ++i) {
        vector[i] = arr[i][j][k];
    }
    return vector;
}

DataType *y_at(int Ny, int i, int k, DataType ***arr) {
    auto *vector = new DataType[Ny];
    for (int j = 0; j < Ny; ++j) {
        vector[j] = arr[i][j][k];
    }
    return vector;
}

DataType *z_at(int Nz, int i, int j, DataType ***arr) {
    auto *vector = new DataType[Nz];
    for (int k = 0; k < Nz; ++k) {
        vector[k] = arr[i][j][k];
    }
    return vector;
}
