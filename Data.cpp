//
// Created by 溫晧良 on 2022/9/28.
//

#include "Data.h"
#include "global.h"
#include <cassert>
#include <iostream>

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

Data::Data(const int nx) {
    ndim = 1;
    _Nx = nx;
    _Ny = 1;
    _Nz = 1;
    data = init_array(_Nx + 2 * ghc, _Ny, _Nz);
}

Data::Data(const int nx, const int ny) {
    ndim = 2;
    _Nx = nx;
    _Ny = ny;
    _Nz = 1;
    data = init_array(_Nx + 2 * ghc, _Ny + 2 * ghc, _Nz);
}

Data::Data(const int nx, const int ny, const int nz) {
    ndim = 3;
    _Nx = nx;
    _Ny = ny;
    _Nz = nz;
    data = init_array(_Nx + 2 * ghc, _Ny + 2 * ghc, _Nz + 2 * ghc);
}

int Data::Nx() {
    return _Nx;
}

int Data::Ny() {
    return _Ny;
}

int Data::Nz() {
    return _Nz;
}

MapIndex Data::index_mapping(int i, int j, int k) {
    unsigned I, J, K;

    bool i_checked = (_Nx > 1 and i > -ghc and i < _Nx + ghc + 1) or _Nx == i;
    bool j_checked = (_Ny > 1 and j > -ghc and j < _Ny + ghc + 1) or _Ny == j;
    bool k_checked = (_Nz > 1 and k > -ghc and k < _Nz + ghc + 1) or _Nz == k;

    assert(i_checked);
    assert(j_checked);
    assert(k_checked);

    I = i - 1;
    J = j - 1;
    K = k - 1;

    if (_Nx > 1) {
        I = I + ghc;
    }

    if (_Ny > 1) {
        J = J + ghc;
    }

    if (_Nz > 1) {
        K = K + ghc;
    }

    auto mapped_index = MapIndex{I, J, K};
    return mapped_index;
}