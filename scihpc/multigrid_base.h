//
// Created by 溫晧良 on 2022/10/14.
//

#ifndef SCIHPC_MULTIGRID_BASE_H
#define SCIHPC_MULTIGRID_BASE_H

#include "global.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::Matrix<DataType, Eigen::Dynamic, 1> VectorX;
typedef Eigen::SparseMatrix<DataType> SMatrix;
typedef Eigen::SparseLU<SMatrix, Eigen::COLAMDOrdering<int> > SpSolver;
typedef Eigen::ConjugateGradient<SMatrix, Eigen::Lower | Eigen::Upper> CGSolver;

class multigrid_base {
public:

    multigrid_base(int _nx, int _ny, int _degree, DataType _dx, DataType _dy);

    multigrid_base(int _nx, int _ny, int _nz, int _degree, DataType _dx, DataType _dy, DataType _dz);

    void init_full();

    SMatrix A;
    VectorX rhs, sol, res, buffer;
    SpSolver *solver;
    CGSolver *cg;

    int nx, ny, nz, degree, ndim, n;
    DataType dx, dy, dz;
    bool no_compatibility = false;

    int base_step = 10;

    int of(int i, int j, int k);

    int of(int i, int j);

    void Ax(const VectorX &x);

    void compatibility_condition(VectorX &f);

    void relax(int iter);

    DataType residual();

    void restriction(multigrid_base &coarse);

    void prolongation(multigrid_base &dense);

    void solve(DataType tol);
};

#endif //SCIHPC_MULTIGRID_BASE_H
