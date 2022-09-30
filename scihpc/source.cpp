//
// Created by 溫晧良 on 2022/10/1.
//

#include "source.h"

void convection(scalar_data *f, vector_data *vel, structured_grid *geo, DataType ***s,
                void(*flux)(scalar_data *,vector_data *)) {
    flux(f, vel);
    uccd_find_derivatives(f, geo, vel);
#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                s[i][j][k] = -vel->x.data[i][j][k] * f->fx[i][j][k];
            }
        }
    }
    if (f->ndim > 1) {
#pragma omp parallel for
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                for (int k = 0; k < f->Nz; ++k) {
                    s[i][j][k] -= vel->y.data[i][j][k] * f->fy[i][j][k];
                }
            }
        }
    }
    if (f->ndim > 2) {
#pragma omp parallel for
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                for (int k = 0; k < f->Nz; ++k) {
                    s[i][j][k] -= vel->z.data[i][j][k] * f->fz[i][j][k];
                }
            }
        }
    }
}

void Hamilton_Jacobi(scalar_data *f, vector_data *vel, structured_grid *geo, DataType ***s,
                     void (*flux)(scalar_data *, vector_data *)) {
    flux(f, vel);
    uccd_find_derivatives(f, geo, vel);
#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                s[i][j][k] = -f->fx[i][j][k];
            }
        }
    }
    if (f->ndim > 1) {
#pragma omp parallel for
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                for (int k = 0; k < f->Nz; ++k) {
                    s[i][j][k] -= f->fy[i][j][k];
                }
            }
        }
    }
    if (f->ndim > 2) {
#pragma omp parallel for
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                for (int k = 0; k < f->Nz; ++k) {
                    s[i][j][k] -= f->fz[i][j][k];
                }
            }
        }
    }

}
