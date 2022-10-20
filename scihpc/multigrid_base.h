//
// Created by 溫晧良 on 2022/10/14.
//

#ifndef SCIHPC_MULTIGRID_BASE_H
#define SCIHPC_MULTIGRID_BASE_H

#include "SparseMatrix.h"

class multigrid_base {
public:
    multigrid_base();

    void init3d(int _nx, int _ny, int _nz, int _degree, DataType _dx, DataType _dy, DataType _dz);

    void init2d(int _nx, int _ny, int _degree, DataType _dx, DataType _dy);

    void init_full();

    SparseMatrix *oD;
    DataType *rhs, *sol, *D, *res;
    int nx, ny, nz, degree, ndim, n;
    DataType dx, dy, dz;
    bool no_compatibility = false;

    int of(int i, int j, int k);

    int of(int i, int j);

    void compatibility_condition(DataType *f);

    void relax(int iter);

    DataType residual();

    void restriction(multigrid_base &coarse);

    void prolongation(multigrid_base &dense);

    void solve(DataType tol);
};

#endif //SCIHPC_MULTIGRID_BASE_H
