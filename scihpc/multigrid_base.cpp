//
// Created by 溫晧良 on 2022/10/14.
//

#include "multigrid_base.h"

template<typename T>
multigrid_base<T>::multigrid_base(int _nx, int _ny, int _nz, int _degree, DataType _dx, DataType _dy, DataType _dz) {
    nx = _nx;
    ny = _ny;
    nz = _nz;
    degree = _degree;
    A = SparseMatrix<T>(nx * ny * nz);
    dx = _dx;
    dy = _dy;
    dz = _dz;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                auto c = 0.0;
                if (i > 0) {
                    A.set(1.0 / (dx * dx), of(i, j, k), of(i - 1, j, k));
                    c += 1.0 / (dx * dx);
                }
                if (i < nx - 1) {
                    A.set(1.0 / (dx * dx), of(i, j, k), of(i + 1, j, k));
                    c += 1.0 / (dx * dx);
                }
                if (j > 0) {
                    A.set(1.0 / (dy * dy), of(i, j, k), of(i, j - 1, k));
                    c += 1.0 / (dy * dy);
                }
                if (j < ny - 1) {
                    A.set(1.0 / (dy * dy), of(i, j, k), of(i, j + 1, k));
                    c += 1.0 / (dy * dy);
                }
                if (k > 0) {
                    A.set(1.0 / (dz * dz), of(i, j, k), of(i, j, k - 1));
                    c += 1.0 / (dz * dz);
                }
                if (k < nz - 1) {
                    A.set(1.0 / (dz * dz), of(i, j, k), of(i, j, k + 1));
                    c += 1.0 / (dz * dz);
                }
                A.set(-c, of(i, j, k), of(i, j, k));
            }
        }
    }

}

template<typename T>
multigrid_base<T>::multigrid_base(int _nx, int _ny, int _degree, DataType _dx, DataType _dy) {
    nx = _nx;
    ny = _ny;
    degree = _degree;
    A = SparseMatrix<T>(nx * ny);
    dx = _dx;
    dy = _dy;
    dz = 0.0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            auto c = 0.0;
            if (i > 0) {
                A.set(1.0 / (dx * dx), of(i, j), of(i - 1, j));
                c += 1.0 / (dx * dx);
            }
            if (i < nx - 1) {
                A.set(1.0 / (dx * dx), of(i, j), of(i + 1, j));
                c += 1.0 / (dx * dx);
            }
            if (j > 0) {
                A.set(1.0 / (dy * dy), of(i, j), of(i, j - 1));
                c += 1.0 / (dy * dy);
            }
            if (j < ny - 1) {
                A.set(1.0 / (dy * dy), of(i, j), of(i, j + 1));
                c += 1.0 / (dy * dy);
            }
            A.set(-c, of(i, j), of(i, j));
        }
    }
}

template<typename T>
int multigrid_base<T>::of(int i, int j, int k) {
    return i + j * nx + k * nx * ny + 1;
}

template<typename T>
int multigrid_base<T>::of(int i, int j) {
    return i + j * nx + 1;
}
