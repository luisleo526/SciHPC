//
// Created by 溫晧良 on 2022/10/14.
//

#ifndef SCIHPC_MULTIGRID_BASE_H
#define SCIHPC_MULTIGRID_BASE_H

#include "global.h"
#include "SparseMatrix.h"

template<typename T>
class multigrid_base {
public:
    multigrid_base(int _nx, int _ny, int _nz, int _degree, DataType _dx, DataType _dy, DataType _dz);
    multigrid_base(int _nx, int _ny, int _degree, DataType _dx, DataType _dy);
    int nx, ny, nz, degree;
    DataType dx, dy, dz;
    SparseMatrix<T> A;
    int of(int i, int j, int k);
    int of(int i, int j);
};

#endif //SCIHPC_MULTIGRID_BASE_H
