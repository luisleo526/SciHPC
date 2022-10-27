//
// Created by 溫晧良 on 2022/10/27.
//

#include "simple_bc.h"

void zero_order_extrapolation(DataType ***f, int nx, int ny, int nz, int ghc) {

#pragma omp parallel for default(none) shared(f, nx, ny, nz, ghc) collapse(3)
    for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < nz; ++k) {
            for (int i = 0; i < ghc; ++i) {
                f[i][j][k] = f[ghc][j][k];
                f[nx - 1 - i][j][k] = f[nx - 1 - ghc][j][k];
            }
        }
    }

    if (ny > 1) {
#pragma omp parallel for default(none) shared(f, nx, ny, nz, ghc) collapse(3)
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) {
                for (int j = 0; j < ghc; ++j) {
                    f[i][j][k] = f[i][ghc][k];
                    f[i][ny - 1 - j][k] = f[i][ny - 1 - ghc][k];
                }
            }
        }
    }

    if (nz > 1) {
#pragma omp parallel for default(none) shared(f, nx, ny, nz, ghc) collapse(3)
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < ghc; ++k) {
                    f[i][j][k] = f[i][j][ghc];
                    f[i][j][nz - 1 - k] = f[i][j][nz - 1 - ghc];
                }
            }
        }
    }

}

void first_order_extrapolation(DataType ***f, int nx, int ny, int nz, int ghc) {

#pragma omp parallel for default(none) shared(f, nx, ny, nz, ghc) collapse(3)
    for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < nz; ++k) {
            for (int i = 0; i < ghc; ++i) {
                f[i][j][k] = 2 * f[ghc][j][k] - f[ghc + 1][j][k];
                f[nx - 1 - i][j][k] = 2 * f[nx - 1 - ghc][j][k] - f[nx - 1 - ghc - 1][j][k];
            }
        }
    }

    if (ny > 1) {
#pragma omp parallel for default(none) shared(f, nx, ny, nz, ghc) collapse(3)
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) {
                for (int j = 0; j < ghc; ++j) {
                    f[i][j][k] = 2 * f[i][ghc][k] - f[i][ghc + 1][k];
                    f[i][ny - 1 - j][k] = 2 * f[i][ny - 1 - ghc][k] - f[i][ny - 1 - ghc - 1][k];
                }
            }
        }
    }

    if (nz > 1) {
#pragma omp parallel for default(none) shared(f, nx, ny, nz, ghc) collapse(3)
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < ghc; ++k) {
                    f[i][j][k] = 2 * f[i][j][ghc] - f[i][j][ghc + 1];
                    f[i][j][nz - 1 - k] = 2 * f[i][j][nz - 1 - ghc] - f[i][j][nz - 1 - ghc - 1];
                }
            }
        }
    }
}
