//
// Created by 溫晧良 on 2022/10/11.
//

#include "second_order_solver.h"

second_order_solver::second_order_solver(DataType _dx, DataType _dy, DataType _dz) {
    dx = _dx;
    dy = _dy;
    dz = _dz;
}

void second_order_solver::find_fx(scalar_data *f) const {
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                if (i == 0) {
                    f->fx[i][j][k] =
                            (-1.5 * f->flux[i][j][k] + 2 * f->flux[i + 1][j][k] - 0.5 * f->flux[i + 2][j][k]) / dx;
                    f->fxx[i][j][k] =
                            (2 * f->flux[i][j][k] - 5 * f->flux[i + 1][j][k] + 4 * f->flux[i + 2][j][k] -
                             f->flux[i + 3][j][k]) / (dx * dx);
                } else if (i == f->Nx - 1) {
                    f->fx[i][j][k] =
                            (1.5 * f->flux[i][j][k] - 2 * f->flux[i - 1][j][k] + 0.5 * f->flux[i - 2][j][k]) / dx;
                    f->fxx[i][j][k] =
                            (2 * f->flux[i][j][k] - 5 * f->flux[i - 1][j][k] + 4 * f->flux[i - 2][j][k] -
                             f->flux[i - 3][j][k]) / (dx * dx);
                } else {
                    f->fx[i][j][k] = (f->flux[i + 1][j][k] - f->flux[i - 1][j][k]) / (2 * dx);
                    f->fxx[i][j][k] = (f->flux[i + 1][j][k] - 2 * f->flux[i][j][k] + f->flux[i - 1][j][k]) / (dx * dx);
                }
            }
        }
    }
}

void second_order_solver::find_fy(scalar_data *f) const {
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                if (j == 0) {
                    f->fy[i][j][k] =
                            (-1.5 * f->flux[i][j][k] + 2 * f->flux[i][j + 1][k] - 0.5 * f->flux[i][j + 2][k]) / dy;
                    f->fyy[i][j][k] =
                            (2 * f->flux[i][j][k] - 5 * f->flux[i][j + 1][k] + 4 * f->flux[i][j + 2][k] -
                             f->flux[i][j + 3][k]) / (dy * dy);
                } else if (j == f->Ny - 1) {
                    f->fy[i][j][k] =
                            (1.5 * f->flux[i][j][k] - 2 * f->flux[i][j - 1][k] + 0.5 * f->flux[i][j - 2][k]) / dy;
                    f->fyy[i][j][k] =
                            (2 * f->flux[i][j][k] - 5 * f->flux[i][j - 1][k] + 4 * f->flux[i][j - 2][k] -
                             f->flux[i][j - 3][k]) / (dy * dy);
                } else {
                    f->fy[i][j][k] = (f->flux[i][j + 1][k] - f->flux[i][j - 1][k]) / (2 * dy);
                    f->fyy[i][j][k] = (f->flux[i][j + 1][k] - 2 * f->flux[i][j][k] + f->flux[i][j - 1][k]) / (dy * dy);
                }
            }
        }
    }
}

void second_order_solver::find_fz(scalar_data *f) const {
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                if (k == 0) {
                    f->fz[i][j][k] =
                            (-1.5 * f->flux[i][j][k] + 2 * f->flux[i][j][k + 1] - 0.5 * f->flux[i][j][k + 2]) / dz;
                    f->fzz[i][j][k] =
                            (2 * f->flux[i][j][k] - 5 * f->flux[i][j][k + 1] + 4 * f->flux[i][j][k + 2] -
                             f->flux[i][j][k + 3]) / (dz * dz);
                } else if (k == f->Nz - 1) {
                    f->fz[i][j][k] =
                            (1.5 * f->flux[i][j][k] - 2 * f->flux[i][j][k - 1] + 0.5 * f->flux[i][j][k - 2]) / dz;
                    f->fzz[i][j][k] =
                            (2 * f->flux[i][j][k] - 5 * f->flux[i][j][k - 1] + 4 * f->flux[i][j][k - 2] -
                             f->flux[i][j][k - 3]) / (dz * dz);
                } else {
                    f->fz[i][j][k] = (f->flux[i][j][k + 1] - f->flux[i][j][k - 1]) / (2 * dz);
                    f->fzz[i][j][k] = (f->flux[i][j][k + 1] - 2 * f->flux[i][j][k] + f->flux[i][j][k - 1]) / (dz * dz);
                }
            }
        }
    }
}

void second_order_solver::find_fx(scalar_data *f, vector_data *vel) const {
#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int I = 0; I < f->nx; ++I) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index = f->index_mapping(I + 1, 1, 1);
                auto i = index.i;

                if (vel->x.data[i][j][k] > 0.0) {
                    f->fx[i][j][k] =
                            (1.5 * f->flux[i][j][k] - 2 * f->flux[i - 1][j][k] + 0.5 * f->flux[i - 2][j][k]) / dx;
                } else {
                    f->fx[i][j][k] =
                            (-1.5 * f->flux[i][j][k] + 2 * f->flux[i + 1][j][k] - 0.5 * f->flux[i + 2][j][k]) / dx;
                }
                f->fxx[i][j][k] = (f->flux[i + 1][j][k] - 2 * f->flux[i][j][k] + f->flux[i - 1][j][k]) / (dx * dx);
            }
        }
    }
}

void second_order_solver::find_fy(scalar_data *f, vector_data *vel) const {
#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int J = 0; J < f->ny; ++J) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index = f->index_mapping(1, J + 1, 1);
                auto j = index.j;
                if (vel->y.data[i][j][k] > 0.0) {
                    f->fy[i][j][k] =
                            (1.5 * f->flux[i][j][k] - 2 * f->flux[i][j - 1][k] + 0.5 * f->flux[i][j - 2][k]) / dy;
                } else {
                    f->fy[i][j][k] =
                            (-1.5 * f->flux[i][j][k] + 2 * f->flux[i][j + 1][k] - 0.5 * f->flux[i][j + 2][k]) / dy;
                }
                f->fyy[i][j][k] = (f->flux[i][j + 1][k] - 2 * f->flux[i][j][k] + f->flux[i][j - 1][k]) / (dy * dy);
            }
        }
    }
}

void second_order_solver::find_fz(scalar_data *f, vector_data *vel) const {
#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int K = 0; K < f->nz; ++K) {
                auto index = f->index_mapping(1, 1, K + 1);
                auto k = index.k;
                if (vel->z.data[i][j][k] > 0.0) {
                    f->fz[i][j][k] =
                            (1.5 * f->flux[i][j][k] - 2 * f->flux[i][j][k - 1] + 0.5 * f->flux[i][j][k - 2]) / dz;
                } else {
                    f->fz[i][j][k] =
                            (-1.5 * f->flux[i][j][k] + 2 * f->flux[i][j][k + 1] - 0.5 * f->flux[i][j][k + 2]) / dz;
                }
                f->fzz[i][j][k] = (f->flux[i][j][k + 1] - 2 * f->flux[i][j][k] + f->flux[i][j][k - 1]) / (dz * dz);
            }
        }
    }
}

void second_order_solver::find_derivatives(scalar_data *f) const {
    find_fx(f);
    if (f->ndim > 1){
        find_fy(f);
    }
    if (f->ndim > 2){
        find_fz(f);
    }
}

void second_order_solver::find_derivatives(scalar_data *f, vector_data *vel) const {
    find_fx(f, vel);
    if (f->ndim > 1){
        find_fy(f, vel);
    }
    if (f->ndim > 2){
        find_fz(f, vel);
    }
}




