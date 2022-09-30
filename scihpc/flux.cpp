//
// Created by 溫晧良 on 2022/10/1.
//

#include "flux.h"

void identity_flux(scalar_data *f, vector_data *vel) {
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->data[i][j][k];
            }
        }
    }
}

void burgers_flux(scalar_data *f, vector_data *vel) {
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->data[i][j][k] * f->data[i][j][k] * 0.5;
            }
        }
    }
}
