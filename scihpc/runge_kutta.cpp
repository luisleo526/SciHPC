//
// Created by 溫晧良 on 2022/9/28.
//

#include "runge_kutta.h"

runge_kutta::runge_kutta(int nx, int ny, int nz) {
    s1 = init_array(nx, ny, nz);
    s2 = init_array(nx, ny, nz);
    s3 = init_array(nx, ny, nz);
}

void
runge_kutta::tvd_rk3(wrapper *f, wrapper *vel, void (*flux)(scalar_data *, vector_data *),
                     void (*rhs)(wrapper *, wrapper *, DataType ***, void (*)(scalar_data *, vector_data *))) const {

    (*rhs)(f, vel, s1, flux);
#pragma omp parallel for default(none) shared(f, s1, geo)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->scalar->data[i][j][k] += f->params->dt * s1[i][j][k];
            }
        }
    }
    f->apply_scalar_bc();

    (*rhs)(f, vel, s2, flux);
#pragma omp parallel for default(none) shared(f, s1, s2, geo)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->scalar->data[i][j][k] += f->params->dt * (s2[i][j][k] - 3.0 * s1[i][j][k]) / 4.0;
            }
        }
    }
    f->apply_scalar_bc();

    (*rhs)(f, vel, s3, flux);
#pragma omp parallel for default(none) shared(f, s1, s2, s3, geo)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->scalar->data[i][j][k] += f->params->dt * (-s1[i][j][k] - s2[i][j][k] + 8.0 * s3[i][j][k]) / 12.0;
            }
        }
    }
    f->apply_scalar_bc();
}

void runge_kutta::euler(wrapper *f, wrapper *vel, void (*flux)(scalar_data *, vector_data *),
                        void (*rhs)(wrapper *, wrapper *, DataType ***, void (*)(scalar_data *, vector_data *))) const {
    (*rhs)(f, vel, s1, flux);
#pragma omp parallel for default(none) shared(f, s1, geo)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->scalar->data[i][j][k] += f->params->dt * s1[i][j][k];
            }
        }
    }
    f->apply_scalar_bc();
}
