//
// Created by leo on 10/2/22.
//

#include "uccd_solver.h"

uccd_solver::uccd_solver(int _nx, int _ny, int _nz, DataType _dx, DataType _dy, DataType _dz) {
    x = new uccd_base(_nx, _dx);
    y = new uccd_base(_ny, _dy);
    z = new uccd_base(_nz, _dz);
}

uccd_solver::uccd_solver(scalar_data *f, structured_grid *geo) {
    x = new uccd_base(f->Nx, geo->dx);
    y = new uccd_base(f->Ny, geo->dy);
    z = new uccd_base(f->Nz, geo->dz);
}


void uccd_solver::find_fx(scalar_data *f, vector_data *vel) const {

    DataType c1, c2, c3;
    c3 = -0.06096119008109;
    c2 = 1.99692238016218;
    c1 = -1.93596119008109;

    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {

            x->upwind->s[0] = (-3.5 * f->flux[0][j][k]
                               + 4.0 * f->flux[1][j][k]
                               - 0.5 * f->flux[2][j][k]) / x->h;

            x->upwind->ss[0] = (34.0 / 3.0 * f->flux[0][j][k]
                                - 83.0 / 4.0 * f->flux[1][j][k]
                                + 10.0 * f->flux[2][j][k]
                                - 7.0 / 12.0 * f->flux[3][j][k]) / x->h / x->h;

            x->upwind->s[f->Nx - 1] = -(-3.5 * f->flux[f->Nx - 1][j][k]
                                        + 4.0 * f->flux[f->Nx - 2][j][k]
                                        - 0.5 * f->flux[f->Nx - 3][j][k]) / x->h;
            x->upwind->ss[f->Nx - 1] = (34.0 / 3.0 * f->flux[f->Nx - 1][j][k]
                                        - 83.0 / 4.0 * f->flux[f->Nx - 2][j][k]
                                        + 10.0 * f->flux[f->Nx - 3][j][k]
                                        - 7.0 / 12.0 * f->flux[f->Nx - 4][j][k]) / x->h / x->h;

            x->downwind->s[0] = x->upwind->s[0];
            x->downwind->ss[0] = x->upwind->ss[0];
            x->downwind->s[f->Nx - 1] = x->upwind->s[f->Nx - 1];
            x->downwind->ss[f->Nx - 1] = x->upwind->ss[f->Nx - 1];

            for (int i = 1; i < f->Nx - 1; ++i) {
                x->upwind->s[i] = (c1 * f->flux[i - 1][j][k]
                                   + c2 * f->flux[i][j][k]
                                   + c3 * f->flux[i + 1][j][k]) / x->h;
                x->upwind->ss[i] = (3.0 * f->flux[i - 1][j][k] - 6.0 * f->flux[i][j][k]
                                    + 3.0 * f->flux[i + 1][j][k]) / x->h / x->h;

                x->downwind->s[i] = -(c3 * f->flux[i - 1][j][k]
                                      + c2 * f->flux[i][j][k]
                                      + c1 * f->flux[i + 1][j][k]) / x->h;
                x->downwind->ss[i] = x->upwind->ss[i];
            }

            twin_bks(x->upwind->a, x->upwind->b,
                     x->upwind->aa, x->upwind->bb,
                     x->upwind->s, x->upwind->ss, f->Nx);
            twin_bks(x->downwind->a, x->downwind->b,
                     x->downwind->aa, x->downwind->bb,
                     x->downwind->s, x->downwind->ss, f->Nx);

            for (int i = 0; i < f->Nx; ++i) {
                if (vel->x.data[i][j][k] > 0) {
                    f->fx[i][j][k] = x->upwind->s[i];
                } else {
                    f->fx[i][j][k] = x->downwind->s[i];
                }
                f->fxx[i][j][k] = 0.5 * (x->upwind->ss[i] + x->downwind->ss[i]);
            }

        }
    }
}

void uccd_solver::find_fy(scalar_data *f, vector_data *vel) const {

    DataType c1, c2, c3;
    c3 = -0.06096119008109;
    c2 = 1.99692238016218;
    c1 = -1.93596119008109;

    for (int i = 0; i < f->Nx; ++i) {
        for (int k = 0; k < f->Nz; ++k) {

            y->upwind->s[0] = (-3.5 * f->flux[i][0][k]
                               + 4.0 * f->flux[i][1][k]
                               - 0.5 * f->flux[i][2][k]) / y->h;

            y->upwind->ss[0] = (34.0 / 3.0 * f->flux[i][0][k]
                                - 83.0 / 4.0 * f->flux[i][1][k]
                                + 10.0 * f->flux[i][2][k]
                                - 7.0 / 12.0 * f->flux[i][3][k]) / y->h / y->h;

            y->upwind->s[f->Ny - 1] = -(-3.5 * f->flux[i][f->Ny - 1][k]
                                        + 4.0 * f->flux[i][f->Ny - 2][k]
                                        - 0.5 * f->flux[i][f->Ny - 3][k]) / y->h;
            y->upwind->ss[f->Ny - 1] = (34.0 / 3.0 * f->flux[i][f->Ny - 1][k]
                                        - 83.0 / 4.0 * f->flux[i][f->Ny - 2][k]
                                        + 10.0 * f->flux[i][f->Ny - 3][k]
                                        - 7.0 / 12.0 * f->flux[i][f->Ny - 4][k]) / y->h / y->h;

            y->downwind->s[0] = y->upwind->s[0];
            y->downwind->ss[0] = y->upwind->ss[0];
            y->downwind->s[f->Ny - 1] = y->upwind->s[f->Ny - 1];
            y->downwind->ss[f->Ny - 1] = y->upwind->ss[f->Ny - 1];

            for (int j = 1; j < f->Ny - 1; ++j) {
                y->upwind->s[j] = (c1 * f->flux[i][j - 1][k]
                                   + c2 * f->flux[i][j][k]
                                   + c3 * f->flux[i][j + 1][k]) / y->h;
                y->upwind->ss[j] = (3.0 * f->flux[i][j - 1][k] - 6.0 * f->flux[i][j][k]
                                    + 3.0 * f->flux[i][j + 1][k]) / y->h / y->h;

                y->downwind->s[j] = -(c3 * f->flux[i][j - 1][k]
                                      + c2 * f->flux[i][j][k]
                                      + c1 * f->flux[i][j + 1][k]) / y->h;
                y->downwind->ss[j] = y->upwind->ss[j];
            }

            twin_bks(y->upwind->a, y->upwind->b,
                     y->upwind->aa, y->upwind->bb,
                     y->upwind->s, y->upwind->ss, f->Ny);

            twin_bks(y->downwind->a, y->downwind->b,
                     y->downwind->aa, y->downwind->bb,
                     y->downwind->s, y->downwind->ss, f->Ny);

            for (int j = 0; j < f->Ny; ++j) {
                if (vel->y.data[i][j][k] > 0) {
                    f->fy[i][j][k] = y->upwind->s[j];
                } else {
                    f->fy[i][j][k] = y->downwind->s[j];
                }
                f->fyy[i][j][k] = 0.5 * (y->upwind->ss[j] + y->downwind->ss[j]);
            }
        }
    }
}

void uccd_solver::find_fz(scalar_data *f, vector_data *vel) const {

    DataType c1, c2, c3;
    c3 = -0.06096119008109;
    c2 = 1.99692238016218;
    c1 = -1.93596119008109;

    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {

            z->upwind->s[0] = (-3.5 * f->flux[i][j][0]
                               + 4.0 * f->flux[i][j][1]
                               - 0.5 * f->flux[i][j][2]) / z->h;

            z->upwind->ss[0] = (34.0 / 3.0 * f->flux[i][j][0]
                                - 83.0 / 4.0 * f->flux[i][j][1]
                                + 10.0 * f->flux[i][j][2]
                                - 7.0 / 12.0 * f->flux[i][j][3]) / z->h / z->h;

            z->upwind->s[f->Nz - 1] = -(-3.5 * f->flux[i][j][f->Nz - 1]
                                        + 4.0 * f->flux[i][j][f->Nz - 2]
                                        - 0.5 * f->flux[i][j][f->Nz - 3]) / z->h;
            z->upwind->ss[f->Nz - 1] = (34.0 / 3.0 * f->flux[i][j][f->Nz - 1]
                                        - 83.0 / 4.0 * f->flux[i][j][f->Nz - 2]
                                        + 10.0 * f->flux[i][j][f->Nz - 3]
                                        - 7.0 / 12.0 * f->flux[i][j][f->Nz - 4]) / z->h / z->h;

            z->downwind->s[0] = z->upwind->s[0];
            z->downwind->ss[0] = z->upwind->ss[0];
            z->downwind->s[f->Nz - 1] = z->upwind->s[f->Nz - 1];
            z->downwind->ss[f->Nz - 1] = z->upwind->ss[f->Nz - 1];

            for (int k = 1; k < f->Nz - 1; ++k) {
                z->upwind->s[k] = (c1 * f->flux[i][j][k - 1]
                                   + c2 * f->flux[i][j][k]
                                   + c3 * f->flux[i][j][k + 1]) / z->h;
                z->upwind->ss[k] = (3.0 * f->flux[i][j][k - 1] - 6.0 * f->flux[i][j][k]
                                    + 3.0 * f->flux[i][j][k + 1]) / z->h / z->h;
                z->downwind->s[k] = -(c3 * f->flux[i][j][k - 1]
                                      + c2 * f->flux[i][j][k]
                                      + c1 * f->flux[i][j][k + 1]) / z->h;
                z->downwind->ss[k] = z->upwind->ss[k];
            }

            twin_bks(z->upwind->a, z->upwind->b,
                     z->upwind->aa, z->upwind->bb,
                     z->upwind->s, z->upwind->ss, f->Nz);

            twin_bks(z->downwind->a, z->downwind->b,
                     z->downwind->aa, z->downwind->bb,
                     z->downwind->s, z->downwind->ss, f->Nz);

            for (int k = 0; k < f->Nz; ++k) {
                if (vel->z.data[i][j][k] > 0) {
                    f->fz[i][j][k] = z->upwind->s[k];
                } else {
                    f->fz[i][j][k] = z->downwind->s[k];
                }
                f->fzz[i][j][k] = 0.5 * (z->upwind->ss[k] + z->downwind->ss[k]);
            }
        }
    }

}

void uccd_solver::find_derivatives(scalar_data *f, vector_data *vel) const {

    find_fx(f, vel);
    if (f->ndim > 1) { find_fy(f, vel); }
    if (f->ndim > 2) { find_fz(f, vel); }
}