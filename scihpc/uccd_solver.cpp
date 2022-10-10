//
// Created by leo on 10/2/22.
//

#include "uccd_solver.h"

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

#pragma omp parallel for default(none) shared(f, vel, c1, c2, c3) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {

            DataType su[f->Nx], ssu[f->Nx], sd[f->Nx], ssd[f->Nx];

            su[0] = (-3.5 * f->flux[0][j][k]
                     + 4.0 * f->flux[1][j][k]
                     - 0.5 * f->flux[2][j][k]) / x->h;

            ssu[0] = (34.0 / 3.0 * f->flux[0][j][k]
                      - 83.0 / 4.0 * f->flux[1][j][k]
                      + 10.0 * f->flux[2][j][k]
                      - 7.0 / 12.0 * f->flux[3][j][k]) / x->h / x->h;

            su[f->Nx - 1] = -(-3.5 * f->flux[f->Nx - 1][j][k]
                              + 4.0 * f->flux[f->Nx - 2][j][k]
                              - 0.5 * f->flux[f->Nx - 3][j][k]) / x->h;
            ssu[f->Nx - 1] = (34.0 / 3.0 * f->flux[f->Nx - 1][j][k]
                              - 83.0 / 4.0 * f->flux[f->Nx - 2][j][k]
                              + 10.0 * f->flux[f->Nx - 3][j][k]
                              - 7.0 / 12.0 * f->flux[f->Nx - 4][j][k]) / x->h / x->h;

            sd[0] = su[0];
            ssd[0] = ssu[0];
            sd[f->Nx - 1] = su[f->Nx - 1];
            ssd[f->Nx - 1] = ssu[f->Nx - 1];

            for (int i = 1; i < f->Nx - 1; ++i) {
                su[i] = (c1 * f->flux[i - 1][j][k]
                         + c2 * f->flux[i][j][k]
                         + c3 * f->flux[i + 1][j][k]) / x->h;
                ssu[i] = (3.0 * f->flux[i - 1][j][k] - 6.0 * f->flux[i][j][k]
                          + 3.0 * f->flux[i + 1][j][k]) / x->h / x->h;

                sd[i] = -(c3 * f->flux[i - 1][j][k]
                          + c2 * f->flux[i][j][k]
                          + c1 * f->flux[i + 1][j][k]) / x->h;
                ssd[i] = ssu[i];
            }

            twin_bks(x->upwind->a, x->upwind->b,
                     x->upwind->aa, x->upwind->bb,
                     su, ssu, f->Nx);
            twin_bks(x->downwind->a, x->downwind->b,
                     x->downwind->aa, x->downwind->bb,
                     sd, ssd, f->Nx);

            for (int i = 0; i < f->Nx; ++i) {
                if (vel->x.data[i][j][k] > 0) {
                    f->fx[i][j][k] = su[i];
                } else {
                    f->fx[i][j][k] = sd[i];
                }
                f->fxx[i][j][k] = 0.5 * (ssu[i] + ssd[i]);
            }
        }
    }
}

void uccd_solver::find_fy(scalar_data *f, vector_data *vel) const {

    DataType c1, c2, c3;
    c3 = -0.06096119008109;
    c2 = 1.99692238016218;
    c1 = -1.93596119008109;

#pragma omp parallel for default(none) shared(f, vel, c1, c2, c3) collapse(2)
    for (int i = 0; i < f->Nx; ++i) {
        for (int k = 0; k < f->Nz; ++k) {

            DataType su[f->Ny], ssu[f->Ny], sd[f->Ny], ssd[f->Ny];

            su[0] = (-3.5 * f->flux[i][0][k]
                     + 4.0 * f->flux[i][1][k]
                     - 0.5 * f->flux[i][2][k]) / y->h;

            ssu[0] = (34.0 / 3.0 * f->flux[i][0][k]
                      - 83.0 / 4.0 * f->flux[i][1][k]
                      + 10.0 * f->flux[i][2][k]
                      - 7.0 / 12.0 * f->flux[i][3][k]) / y->h / y->h;

            su[f->Ny - 1] = -(-3.5 * f->flux[i][f->Ny - 1][k]
                              + 4.0 * f->flux[i][f->Ny - 2][k]
                              - 0.5 * f->flux[i][f->Ny - 3][k]) / y->h;
            ssu[f->Ny - 1] = (34.0 / 3.0 * f->flux[i][f->Ny - 1][k]
                              - 83.0 / 4.0 * f->flux[i][f->Ny - 2][k]
                              + 10.0 * f->flux[i][f->Ny - 3][k]
                              - 7.0 / 12.0 * f->flux[i][f->Ny - 4][k]) / y->h / y->h;

            sd[0] = su[0];
            ssd[0] = ssu[0];
            sd[f->Ny - 1] = su[f->Ny - 1];
            ssd[f->Ny - 1] = ssu[f->Ny - 1];

            for (int j = 1; j < f->Ny - 1; ++j) {
                su[j] = (c1 * f->flux[i][j - 1][k]
                         + c2 * f->flux[i][j][k]
                         + c3 * f->flux[i][j + 1][k]) / y->h;
                ssu[j] = (3.0 * f->flux[i][j - 1][k] - 6.0 * f->flux[i][j][k]
                          + 3.0 * f->flux[i][j + 1][k]) / y->h / y->h;

                sd[j] = -(c3 * f->flux[i][j - 1][k]
                          + c2 * f->flux[i][j][k]
                          + c1 * f->flux[i][j + 1][k]) / y->h;
                ssd[j] = ssu[j];
            }

            twin_bks(y->upwind->a, y->upwind->b,
                     y->upwind->aa, y->upwind->bb,
                     su, ssu, f->Ny);

            twin_bks(y->downwind->a, y->downwind->b,
                     y->downwind->aa, y->downwind->bb,
                     sd, ssd, f->Ny);

            for (int j = 0; j < f->Ny; ++j) {
                if (vel->y.data[i][j][k] > 0) {
                    f->fy[i][j][k] = su[j];
                } else {
                    f->fy[i][j][k] = sd[j];
                }
                f->fyy[i][j][k] = 0.5 * (ssu[j] + ssd[j]);
            }
        }
    }
}

void uccd_solver::find_fz(scalar_data *f, vector_data *vel) const {

    DataType c1, c2, c3;
    c3 = -0.06096119008109;
    c2 = 1.99692238016218;
    c1 = -1.93596119008109;

#pragma omp parallel for default(none) shared(f, vel, c1, c2, c3) collapse(2)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {

            DataType su[f->Nz], ssu[f->Nz], sd[f->Nz], ssd[f->Nz];

            su[0] = (-3.5 * f->flux[i][j][0]
                     + 4.0 * f->flux[i][j][1]
                     - 0.5 * f->flux[i][j][2]) / z->h;

            ssu[0] = (34.0 / 3.0 * f->flux[i][j][0]
                      - 83.0 / 4.0 * f->flux[i][j][1]
                      + 10.0 * f->flux[i][j][2]
                      - 7.0 / 12.0 * f->flux[i][j][3]) / z->h / z->h;

            su[f->Nz - 1] = -(-3.5 * f->flux[i][j][f->Nz - 1]
                              + 4.0 * f->flux[i][j][f->Nz - 2]
                              - 0.5 * f->flux[i][j][f->Nz - 3]) / z->h;
            ssu[f->Nz - 1] = (34.0 / 3.0 * f->flux[i][j][f->Nz - 1]
                              - 83.0 / 4.0 * f->flux[i][j][f->Nz - 2]
                              + 10.0 * f->flux[i][j][f->Nz - 3]
                              - 7.0 / 12.0 * f->flux[i][j][f->Nz - 4]) / z->h / z->h;

            sd[0] = su[0];
            ssd[0] = ssu[0];
            sd[f->Nz - 1] = su[f->Nz - 1];
            ssd[f->Nz - 1] = ssu[f->Nz - 1];

            for (int k = 1; k < f->Nz - 1; ++k) {
                su[k] = (c1 * f->flux[i][j][k - 1]
                         + c2 * f->flux[i][j][k]
                         + c3 * f->flux[i][j][k + 1]) / z->h;
                ssu[k] = (3.0 * f->flux[i][j][k - 1] - 6.0 * f->flux[i][j][k]
                          + 3.0 * f->flux[i][j][k + 1]) / z->h / z->h;
                sd[k] = -(c3 * f->flux[i][j][k - 1]
                          + c2 * f->flux[i][j][k]
                          + c1 * f->flux[i][j][k + 1]) / z->h;
                ssd[k] = ssu[k];
            }

            twin_bks(z->upwind->a, z->upwind->b,
                     z->upwind->aa, z->upwind->bb,
                     su, ssu, f->Nz);

            twin_bks(z->downwind->a, z->downwind->b,
                     z->downwind->aa, z->downwind->bb,
                     sd, ssd, f->Nz);

            for (int k = 0; k < f->Nz; ++k) {
                if (vel->z.data[i][j][k] > 0) {
                    f->fz[i][j][k] = su[k];
                } else {
                    f->fz[i][j][k] = sd[k];
                }
                f->fzz[i][j][k] = 0.5 * (ssu[k] + ssd[k]);
            }
        }
    }

}

void uccd_solver::find_derivatives(scalar_data *f, vector_data *vel) const {

    find_fx(f, vel);
    if (f->ndim > 1) { find_fy(f, vel); }
    if (f->ndim > 2) { find_fz(f, vel); }
}