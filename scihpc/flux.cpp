//
// Created by 溫晧良 on 2022/10/1.
//

#include "flux.h"

void identity_flux(scalar_data *f, vector_data *vel) {

#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->data[i][j][k];
            }
        }
    }
}

void identity_flux(scalar_data *f) {

#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->data[i][j][k];
            }
        }
    }

}

void identity_flux(vector_data *vel) {

#pragma omp parallel for default(none) shared(vel) collapse(3)
    for (int i = 0; i < vel->x.Nx; ++i) {
        for (int j = 0; j < vel->x.Ny; ++j) {
            for (int k = 0; k < vel->x.Nz; ++k) {
                vel->x.flux[i][j][k] = vel->x.data[i][j][k];
                vel->y.flux[i][j][k] = vel->y.data[i][j][k];
                vel->z.flux[i][j][k] = vel->z.data[i][j][k];
            }
        }
    }
}

void burgers_flux(scalar_data *f, vector_data *vel) {

#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->data[i][j][k] * f->data[i][j][k] * 0.5;
            }
        }
    }
}
