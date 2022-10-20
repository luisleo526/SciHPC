//
// Created by 溫晧良 on 2022/10/14.
//

#ifndef SCIHPC_MULTIGRID_BASE_H
#define SCIHPC_MULTIGRID_BASE_H

#include "global.h"

class multigrid_base {
public:

    multigrid_base(int _nx, int _ny, int _degree, DataType _dx, DataType _dy);

    multigrid_base(int _nx, int _ny, int _nz, int _degree, DataType _dx, DataType _dy, DataType _dz);

    void init_full() const;

    DataType *rhs, *sol, *res, *buffer;
    DataType *cc, *cr, *cl, *cu, *cd, *cb, *cf;
    int nx, ny, nz, degree, ndim, n;
    DataType dx, dy, dz;
    bool no_compatibility = false;

    int base_step = 20;

    int of(int i, int j, int k);

    int of(int i, int j);

    void Ax(const DataType *x);

    void compatibility_condition(DataType *f);

    void relax(int iter);

    DataType residual();

    void restriction(multigrid_base &coarse);

    void prolongation(multigrid_base &dense);

    void solve(DataType tol);
};

#endif //SCIHPC_MULTIGRID_BASE_H
