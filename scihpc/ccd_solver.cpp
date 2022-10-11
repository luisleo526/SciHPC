//
// Created by leo on 10/2/22.
//

#include "ccd_solver.h"
#include <iostream>

ccd_solver::ccd_solver(scalar_data *f, structured_grid *geo) {
    x = new ccd_base(f->Nx, geo->dx);
    y = new ccd_base(f->Ny, geo->dy);
    z = new ccd_base(f->Nz, geo->dz);
}

void ccd_solver::find_fx(scalar_data *f) const {

#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {

            DataType s[f->Nx], ss[f->Nx];

            s[0] = (-3.5 * f->flux[0][j][k] + 4.0 * f->flux[1][j][k] - 0.5 * f->flux[2][j][k]) / x->h;

            ss[0] = (34.0 / 3.0 * f->flux[0][j][k] - 83.0 / 4.0 * f->flux[1][j][k]
                        + 10.0 * f->flux[2][j][k] - 7.0 / 12.0 * f->flux[3][j][k]) / x->h / x->h;

            s[f->Nx - 1] = -(-3.5 * f->flux[f->Nx - 1][j][k] + 4.0 * f->flux[f->Nx - 2][j][k]
                                - 0.5 * f->flux[f->Nx - 3][j][k]) / x->h;

            ss[f->Nx - 1] = (34.0 / 3.0 * f->flux[f->Nx - 1][j][k] - 83.0 / 4.0 * f->flux[f->Nx - 2][j][k]
                                + 10.0 * f->flux[f->Nx - 3][j][k] - 7.0 / 12.0 * f->flux[f->Nx - 4][j][k])
                               / x->h / x->h;

            for (int i = 1; i < f->Nx - 1; ++i) {
                s[i] = 15.0 / 16.0 * (f->flux[i + 1][j][k] - f->flux[i - 1][j][k]) / x->h;
                ss[i] = (3.0 * f->flux[i - 1][j][k] - 6.0 * f->flux[i][j][k] + 3.0 * f->flux[i + 1][j][k])
                           / x->h / x->h;
            }

            twin_bks(x->a, x->b, x->aa, x->bb, s, ss, f->Nx);

            for (int i = 0; i < f->Nx; ++i) {
                f->fx[i][j][k] = s[i];
                f->fxx[i][j][k] = ss[i];
            }
        }
    }

}

void ccd_solver::find_fy(scalar_data *f) const {

#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int i = 0; i < f->Nx; ++i) {
        for (int k = 0; k < f->Nz; ++k) {

            DataType s[f->Ny], ss[f->Ny];

            s[0] = (-3.5 * f->flux[i][0][k] + 4.0 * f->flux[i][1][k] - 0.5 * f->flux[i][2][k]) / y->h;

            ss[0] = (34.0 / 3.0 * f->flux[i][0][k] - 83.0 / 4.0 * f->flux[i][1][k]
                        + 10.0 * f->flux[i][2][k] - 7.0 / 12.0 * f->flux[i][3][k]) / y->h / y->h;

            s[f->Ny - 1] = -(-3.5 * f->flux[i][f->Ny - 1][k] + 4.0 * f->flux[i][f->Ny - 2][k]
                                - 0.5 * f->flux[i][f->Ny - 3][k]) / y->h;

            ss[f->Ny - 1] = (34.0 / 3.0 * f->flux[i][f->Ny - 1][k] - 83.0 / 4.0 * f->flux[i][f->Ny - 2][k]
                                + 10.0 * f->flux[i][f->Ny - 3][k] - 7.0 / 12.0 * f->flux[i][f->Ny - 4][k])
                               / y->h / y->h;

            for (int j = 1; j < f->Ny - 1; ++j) {
                s[j] = 15.0 / 16.0 * (f->flux[i][j + 1][k] - f->flux[i][j - 1][k]) / y->h;
                ss[j] = (3.0 * f->flux[i][j - 1][k] - 6.0 * f->flux[i][j][k] + 3.0 * f->flux[i][j + 1][k])
                           / y->h / y->h;
            }

            twin_bks(y->a, y->b, y->aa, y->bb, s, ss, f->Ny);

            for (int j = 0; j < f->Ny; ++j) {
                f->fy[i][j][k] = s[j];
                f->fyy[i][j][k] = ss[j];
            }
        }
    }
}

void ccd_solver::find_fz(scalar_data *f) const {

#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {

            DataType s[f->Nz], ss[f->Nz];

            s[0] = (-3.5 * f->flux[i][j][0] + 4.0 * f->flux[i][j][1] - 0.5 * f->flux[i][j][2]) / z->h;

            ss[0] = (34.0 / 3.0 * f->flux[i][j][0] - 83.0 / 4.0 * f->flux[i][j][1]
                        + 10.0 * f->flux[i][j][2] - 7.0 / 12.0 * f->flux[i][j][3]) / z->h / z->h;

            s[f->Nz - 1] = -(-3.5 * f->flux[i][j][f->Nz - 1] + 4.0 * f->flux[i][j][f->Nz - 2]
                                - 0.5 * f->flux[i][j][f->Nz - 3]) / z->h;

            ss[f->Nz - 1] = (34.0 / 3.0 * f->flux[i][j][f->Nz - 1] - 83.0 / 4.0 * f->flux[i][j][f->Nz - 2]
                                + 10.0 * f->flux[i][j][f->Nz - 3] - 7.0 / 12.0 * f->flux[i][j][f->Nz - 4])
                               / z->h / z->h;

            for (int k = 1; k < f->Nz - 1; ++k) {
                s[k] = 15.0 / 16.0 * (f->flux[i][j][k + 1] - f->flux[i][j][k - 1]) / z->h;
                ss[k] = (3.0 * f->flux[i][j][k - 1] - 6.0 * f->flux[i][j][k] + 3.0 * f->flux[i][j][k + 1])
                           / z->h / z->h;
            }

            twin_bks(z->a, z->b, z->aa, z->bb, s, ss, f->Nz);

            for (int k = 0; k < f->Nz; ++k) {
                f->fz[i][j][k] = s[k];
                f->fzz[i][j][k] = ss[k];
            }
        }
    }
}

void ccd_solver::mixed_xy(scalar_data *f) const {

#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int i = 0; i < f->Nx; ++i) {
        for (int k = 0; k < f->Nz; ++k) {

            DataType s[f->Ny], ss[f->Ny];

            // Boundary condition
            for (int j = 0; j < f->ghc; ++j) {
                f->fx[i][j][k] = f->fx[i][f->ghc][k];
                f->fx[i][f->Ny - 1 - j][k] = f->fx[i][f->Ny - 1 - f->ghc][k];
            }

            s[0] = (-3.5 * f->fx[i][0][k] + 4.0 * f->fx[i][1][k] - 0.5 * f->fx[i][2][k]) / y->h;

            ss[0] = (34.0 / 3.0 * f->fx[i][0][k] - 83.0 / 4.0 * f->fx[i][1][k]
                        + 10.0 * f->fx[i][2][k] - 7.0 / 12.0 * f->fx[i][3][k]) / y->h / y->h;

            s[f->Ny - 1] = -(-3.5 * f->fx[i][f->Ny - 1][k] + 4.0 * f->fx[i][f->Ny - 2][k]
                                - 0.5 * f->fx[i][f->Ny - 3][k]) / y->h;

            ss[f->Ny - 1] = (34.0 / 3.0 * f->fx[i][f->Ny - 1][k] - 83.0 / 4.0 * f->fx[i][f->Ny - 2][k]
                                + 10.0 * f->fx[i][f->Ny - 3][k] - 7.0 / 12.0 * f->fx[i][f->Ny - 4][k])
                               / y->h / y->h;

            for (int j = 1; j < f->Ny - 1; ++j) {
                s[j] = 15.0 / 16.0 * (f->fx[i][j + 1][k] - f->fx[i][j - 1][k]) / y->h;
                ss[j] = (3.0 * f->fx[i][j - 1][k] - 6.0 * f->fx[i][j][k] + 3.0 * f->fx[i][j + 1][k])
                           / y->h / y->h;
            }

            twin_bks(y->a, y->b, y->aa, y->bb, s, ss, f->Ny);

            for (int j = 0; j < f->Ny; ++j) {
                f->fxy[i][j][k] = s[j];
            }
        }
    }
}

void ccd_solver::mixed_yz(scalar_data *f) const {

#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int i = 0; i < f->Nx; ++i) {

            DataType s[f->Nz], ss[f->Nz];

            // Boundary condition
            for (int k = 0; k < f->ghc; ++k) {
                f->fy[i][j][k] = f->fy[i][j][f->ghc];
                f->fy[i][j][f->Nz - 1 - k] = f->fy[i][j][f->Nz - 1 - f->ghc];
            }

            s[0] = (-3.5 * f->fy[i][j][0] + 4.0 * f->fy[i][j][1] - 0.5 * f->fy[i][j][2]) / z->h;

            ss[0] = (34.0 / 3.0 * f->fy[i][j][0] - 83.0 / 4.0 * f->fy[i][j][1]
                        + 10.0 * f->fy[i][j][2] - 7.0 / 12.0 * f->fy[i][j][3]) / z->h / z->h;

            s[f->Nz - 1] = -(-3.5 * f->fy[i][j][f->Nz - 1] + 4.0 * f->fy[i][j][f->Nz - 2]
                                - 0.5 * f->fy[i][j][f->Nz - 3]) / z->h;

            ss[f->Nz - 1] = (34.0 / 3.0 * f->fy[i][j][f->Nz - 1] - 83.0 / 4.0 * f->fy[i][j][f->Nz - 2]
                                + 10.0 * f->fy[i][j][f->Nz - 3] - 7.0 / 12.0 * f->fy[i][j][f->Nz - 4])
                               / z->h / z->h;

            for (int k = 1; k < f->Nz - 1; ++k) {
                s[k] = 15.0 / 16.0 * (f->fy[i][j][k + 1] - f->fy[i][j][k - 1]) / z->h;
                ss[k] = (3.0 * f->fy[i][j][k - 1] - 6.0 * f->fy[i][j][k] + 3.0 * f->fy[i][j][k + 1])
                           / z->h / z->h;
            }

            twin_bks(z->a, z->b, z->aa, z->bb, s, ss, f->Nz);

            for (int k = 0; k < f->Nz; ++k) {
                f->fyz[i][j][k] = s[k];
            }
        }
    }
}

void ccd_solver::mixed_zx(scalar_data *f) const {

#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int k = 0; k < f->Nz; ++k) {
        for (int j = 0; j < f->Ny; ++j) {

            DataType s[f->Nx], ss[f->Nx];

            // Boundary condition
            for (int i = 0; i < f->ghc; ++i) {
                f->fz[i][j][k] = f->fz[f->ghc][j][k];
                f->fz[f->Nx - 1 - i][j][k] = f->fz[f->Nx - 1 - f->ghc][j][k];
            }

            s[0] = (-3.5 * f->fz[0][j][k] + 4.0 * f->fz[1][j][k] - 0.5 * f->fz[2][j][k]) / x->h;

            ss[0] = (34.0 / 3.0 * f->fz[0][j][k] - 83.0 / 4.0 * f->fz[1][j][k]
                        + 10.0 * f->fz[2][j][k] - 7.0 / 12.0 * f->fz[3][j][k]) / x->h / x->h;

            s[f->Nx - 1] = -(-3.5 * f->fz[f->Nx - 1][j][k] + 4.0 * f->fz[f->Nx - 2][j][k]
                                - 0.5 * f->fz[f->Nx - 3][j][k]) / x->h;

            ss[f->Nx - 1] = (34.0 / 3.0 * f->fz[f->Nx - 1][j][k] - 83.0 / 4.0 * f->fz[f->Nx - 2][j][k]
                                + 10.0 * f->fz[f->Nx - 3][j][k] - 7.0 / 12.0 * f->fz[f->Nx - 4][j][k])
                               / x->h / x->h;

            for (int i = 1; i < f->Nx - 1; ++i) {
                s[i] = 15.0 / 16.0 * (f->fz[i + 1][j][k] - f->fz[i - 1][j][k]) / x->h;
                ss[i] = (3.0 * f->fz[i - 1][j][k] - 6.0 * f->fz[i][j][k] + 3.0 * f->fz[i + 1][j][k])
                           / x->h / x->h;
            }

            twin_bks(x->a, x->b, x->aa, x->bb, s, ss, f->Nx);

            for (int i = 0; i < f->Nx; ++i) {
                f->fzx[i][j][k] = s[i];
            }
        }
    }
}

void ccd_solver::find_derivatives(scalar_data *f) const {
    find_fx(f);
    if (f->ndim > 1) { find_fy(f); }
    if (f->ndim > 2) { find_fz(f); }
}

void ccd_solver::find_derivatives(vector_data *f) const {
    find_derivatives(&f->x);
    if (f->x.ndim > 1) { find_derivatives(&f->y); }
    if (f->x.ndim > 2) { find_derivatives(&f->z); }
}

void ccd_solver::find_derivatives_all(scalar_data *f) const {
    find_derivatives(f);
    if (f->ndim > 1) {
        mixed_xy(f);
    }
    if (f->ndim > 2) {
        mixed_yz(f);
        mixed_zx(f);
    }
}

void ccd_solver::find_derivatives_all(vector_data *f) const {
    find_derivatives_all(&f->x);
    if (f->x.ndim > 1) { find_derivatives_all(&f->y); }
    if (f->x.ndim > 2) { find_derivatives_all(&f->z); }
}

