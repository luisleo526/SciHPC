//
// Created by leo on 10/2/22.
//

#include "ccd_solver.h"
#include <iostream>
ccd_solver::ccd_solver(int _nx, int _ny, int _nz, DataType _dx, DataType _dy, DataType _dz) {
    x = new ccd_base(_nx, _dx);
    y = new ccd_base(_ny, _dy);
    z = new ccd_base(_nz, _dz);
}

ccd_solver::ccd_solver(scalar_data *f, structured_grid *geo) {
    x = new ccd_base(f->Nx, geo->dx);
    y = new ccd_base(f->Ny, geo->dy);
    z = new ccd_base(f->Nz, geo->dz);
}

void ccd_solver::find_fx(scalar_data *f) const {

    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            x->s[0] = (-3.5 * f->flux[0][j][k] + 4.0 * f->flux[1][j][k] - 0.5 * f->flux[2][j][k]) / x->h;

            x->ss[0] = (34.0 / 3.0 * f->flux[0][j][k] - 83.0 / 4.0 * f->flux[1][j][k]
                        + 10.0 * f->flux[2][j][k] - 7.0 / 12.0 * f->flux[3][j][k]) / x->h / x->h;

            x->s[f->Nx - 1] = -(-3.5 * f->flux[f->Nx - 1][j][k] + 4.0 * f->flux[f->Nx - 2][j][k]
                                - 0.5 * f->flux[f->Nx - 3][j][k]) / x->h;

            x->ss[f->Nx - 1] = (34.0 / 3.0 * f->flux[f->Nx - 1][j][k] - 83.0 / 4.0 * f->flux[f->Nx - 2][j][k]
                                + 10.0 * f->flux[f->Nx - 3][j][k] - 7.0 / 12.0 * f->flux[f->Nx - 4][j][k])
                               / x->h / x->h;

            for (int i = 1; i < f->Nx - 1; ++i) {
                x->s[i] = 15.0 / 16.0 * (f->flux[i + 1][j][k] - f->flux[i - 1][j][k]) / x->h;
                x->ss[i] = (3.0 * f->flux[i - 1][j][k] - 6.0 * f->flux[i][j][k] + 3.0 * f->flux[i + 1][j][k])
                           / x->h / x->h;
            }

            twin_bks(x->a, x->b, x->aa, x->bb, x->s, x->ss, f->Nx);

            for (int i = 0; i < f->Nx; ++i) {
                f->fx[i][j][k] = x->s[i];
                f->fxx[i][j][k] = x->ss[i];
            }
        }
    }

}

void ccd_solver::find_fy(scalar_data *f) const {

    for (int i = 0; i < f->Nx; ++i) {
        for (int k = 0; k < f->Nz; ++k) {
            y->s[0] = (-3.5 * f->flux[i][0][k] + 4.0 * f->flux[i][1][k] - 0.5 * f->flux[i][2][k]) / y->h;

            y->ss[0] = (34.0 / 3.0 * f->flux[i][0][k] - 83.0 / 4.0 * f->flux[i][1][k]
                        + 10.0 * f->flux[i][2][k] - 7.0 / 12.0 * f->flux[i][3][k]) / y->h / y->h;

            y->s[f->Ny - 1] = -(-3.5 * f->flux[i][f->Ny - 1][k] + 4.0 * f->flux[i][f->Ny - 2][k]
                                - 0.5 * f->flux[i][f->Ny - 3][k]) / y->h;

            y->ss[f->Ny - 1] = (34.0 / 3.0 * f->flux[i][f->Ny - 1][k] - 83.0 / 4.0 * f->flux[i][f->Ny - 2][k]
                                + 10.0 * f->flux[i][f->Ny - 3][k] - 7.0 / 12.0 * f->flux[i][f->Ny - 4][k])
                               / y->h / y->h;

            for (int j = 1; j < f->Ny - 1; ++j) {
                y->s[j] = 15.0 / 16.0 * (f->flux[i][j + 1][k] - f->flux[i][j - 1][k]) / y->h;
                y->ss[j] = (3.0 * f->flux[i][j - 1][k] - 6.0 * f->flux[i][j][k] + 3.0 * f->flux[i][j + 1][k])
                           / y->h / y->h;
            }

            twin_bks(y->a, y->b, y->aa, y->bb, y->s, y->ss, f->Ny);

            for (int j = 0; j < f->Ny; ++j) {
                f->fy[i][j][k] = y->s[j];
                f->fyy[i][j][k] = y->ss[j];
            }
        }
    }
}

void ccd_solver::find_fz(scalar_data *f) const {

    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            z->s[0] = (-3.5 * f->flux[i][j][0] + 4.0 * f->flux[i][j][1] - 0.5 * f->flux[i][j][2]) / z->h;

            z->ss[0] = (34.0 / 3.0 * f->flux[i][j][0] - 83.0 / 4.0 * f->flux[i][j][1]
                        + 10.0 * f->flux[i][j][2] - 7.0 / 12.0 * f->flux[i][j][3]) / z->h / z->h;

            z->s[f->Nz - 1] = -(-3.5 * f->flux[i][j][f->Nz - 1] + 4.0 * f->flux[i][j][f->Nz - 2]
                                - 0.5 * f->flux[i][j][f->Nz - 3]) / z->h;

            z->ss[f->Nz - 1] = (34.0 / 3.0 * f->flux[i][j][f->Nz - 1] - 83.0 / 4.0 * f->flux[i][j][f->Nz - 2]
                                + 10.0 * f->flux[i][j][f->Nz - 3] - 7.0 / 12.0 * f->flux[i][j][f->Nz - 4])
                               / z->h / z->h;

            for (int k = 1; k < f->Nz - 1; ++k) {
                z->s[k] = 15.0 / 16.0 * (f->flux[i][j][k + 1] - f->flux[i][j][k - 1]) / z->h;
                z->ss[k] = (3.0 * f->flux[i][j][k - 1] - 6.0 * f->flux[i][j][k] + 3.0 * f->flux[i][j][k + 1])
                           / z->h / z->h;
            }

            twin_bks(z->a, z->b, z->aa, z->bb, z->s, z->ss, f->Nz);

            for (int k = 0; k < f->Nz; ++k) {
                f->fz[i][j][k] = z->s[k];
                f->fzz[i][j][k] = z->ss[k];
            }
        }
    }
}

void ccd_solver::find_derivatives(scalar_data *f) const {
    find_fx(f);
    if (f->ndim > 1) { find_fy(f); }
    if (f->ndim > 2) { find_fz(f); }
}
