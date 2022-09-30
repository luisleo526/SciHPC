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
runge_kutta::tvd_rk3(DataType dt, scalar_data *f, vector_data *vel, structured_grid *geo,
                     void(*flux)(scalar_data *, vector_data *pData),
                     void (*bc)(scalar_data *),
                     void (*rhs)(scalar_data *, vector_data *, structured_grid *, DataType ***,
                                 void (*)(scalar_data *, vector_data *))) {
    (*rhs)(f, vel, geo, s1, flux);
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->data[i][j][k] += dt * s1[i][j][k];
            }
        }
    }
    (*bc)(f);

    (*rhs)(f, vel, geo, s2, flux);
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->data[i][j][k] += dt * (s2[i][j][k] - 3.0 * s1[i][j][k]) / 4.0;
            }
        }
    }
    (*bc)(f);

    (*rhs)(f, vel, geo, s3, flux);
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->data[i][j][k] += dt * (-s1[i][j][k] - s2[i][j][k] + 8.0 * s3[i][j][k]) / 12.0;
            }
        }
    }
    (*bc)(f);
}
