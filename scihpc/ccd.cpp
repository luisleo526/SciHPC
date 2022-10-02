//
// Created by 溫晧良 on 2022/9/30.
//

#include "ccd.h"

void ccd_coefficient_boundary_condition(DataType ***coeff, const int n, const DataType dx) {

    coeff[0][0][0] = 0.0;
    coeff[1][0][0] = 0.0;
    coeff[2][0][0] = 0.0;
    coeff[3][0][0] = 0.0;

    coeff[0][1][0] = 1.0;
    coeff[1][1][0] = 0.0;
    coeff[2][1][0] = 0.0;
    coeff[3][1][0] = 1.0;

    coeff[0][2][0] = 2.0;
    coeff[1][2][0] = -dx;
    coeff[2][2][0] = -2.5 / dx;
    coeff[3][2][0] = 8.5;

    coeff[0][0][n - 1] = 2.0;
    coeff[1][0][n - 1] = dx;
    coeff[2][0][n - 1] = 2.5 / dx;
    coeff[3][0][n - 1] = 8.5;

    coeff[0][1][n - 1] = 1.0;
    coeff[1][1][n - 1] = 0.0;
    coeff[2][1][n - 1] = 0.0;
    coeff[3][1][n - 1] = 1.0;

    coeff[0][2][n - 1] = 0.0;
    coeff[1][2][n - 1] = 0.0;
    coeff[2][2][n - 1] = 0.0;
    coeff[3][2][n - 1] = 0.0;

}

DataType ***ccd_coefficient_matrix(const int n, const DataType h) {

    auto coeff = new DataType **[4]; // [a, b, aa, bb]
    for (int i = 0; i < 4; ++i) {
        coeff[i] = new DataType *[3]; // [i-1, i, i+1]
        for (int j = 0; j < 3; ++j) {
            coeff[i][j] = new DataType[n];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        coeff[0][0][i] = 7.0 / 16.0;
        coeff[0][1][i] = 1.0;
        coeff[0][2][i] = 7.0 / 16.0;

        coeff[1][0][i] = h / 16.0;
        coeff[1][1][i] = 0.0;
        coeff[1][2][i] = -h / 16.0;

        coeff[2][0][i] = -9.0 / 8.0 / h;
        coeff[2][1][i] = 0.0;
        coeff[2][2][i] = 9.0 / 8.0 / h;

        coeff[3][0][i] = -1.0 / 8.0;
        coeff[3][1][i] = 1.0;
        coeff[3][2][i] = -1.0 / 8.0;
    }

    ccd_coefficient_boundary_condition(coeff, n, h);
    twin_dec(coeff[0], coeff[1], coeff[2], coeff[3], n);

    return coeff;
}

void ccd_find_fx(scalar_data *f, structured_grid *geo) {

#pragma omp parallel for
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k <= f->Nz; ++k) {

            auto coeff = ccd_coefficient_matrix(f->Nx, geo->dx);

            // Find src
            DataType *s, *ss;
            s = new DataType[f->Nx];
            ss = new DataType[f->Nx];

            // Src boundary condition
            s[0] = (-3.5 * f->flux[0][j][k] + 4.0 * f->flux[1][j][k] - 0.5 * f->flux[2][j][k]) / geo->dx;

            ss[0] = (34.0 / 3.0 * f->flux[0][j][k] - 83.0 / 4.0 * f->flux[1][j][k]
                     + 10.0 * f->flux[2][j][k] - 7.0 / 12.0 * f->flux[3][j][k]) / geo->dx / geo->dx;

            s[f->Nx - 1] = -(-3.5 * f->flux[f->Nx - 1][j][k] + 4.0 * f->flux[f->Nx - 2][j][k]
                             - 0.5 * f->flux[f->Nx - 3][j][k]) / geo->dx;

            ss[f->Nx - 1] = (34.0 / 3.0 * f->flux[f->Nx - 1][j][k] - 83.0 / 4.0 * f->flux[f->Nx - 2][j][k]
                             + 10.0 * f->flux[f->Nx - 3][j][k] - 7.0 / 12.0 * f->flux[f->Nx - 4][j][k])
                            / geo->dx / geo->dx;

            for (int i = 1; i < f->Nx - 1; ++i) {
                s[i] = 15.0 / 16.0 * (f->flux[i + 1][j][k] - f->flux[i - 1][j][k]) / geo->dx;
                ss[i] = (3.0 * f->flux[i - 1][j][k] - 6.0 * f->flux[i][j][k] + 3.0 * f->flux[i + 1][j][k])
                        / geo->dx / geo->dx;
            }

            twin_bks(coeff[0], coeff[1], coeff[2], coeff[3], s, ss, f->Nx);

            for (int i = 0; i < f->Nx; ++i) {
                f->fx[i][j][k] = s[i];
                f->fxx[i][j][k] = ss[i];
            }

            delete3d(coeff, 4, 3);
            delete[] s;
            delete[] ss;
        }
    }
}


void ccd_find_fy(scalar_data *f, structured_grid *geo) {


}


void ccd_find_derivatives(scalar_data *f, structured_grid *geo) {
    ccd_find_fx(f, geo);
    if (f->ndim > 1) {
        ccd_find_fy(f, geo);
    }
    if (f->ndim > 2) {
        ccd_find_fz(f, geo);
    }
}

DataType ****uccd_coefficient_matrix(const int n, const DataType h) {

    auto coeff = new DataType ***[2]; // upwind, downwind
    for (int i = 0; i < 2; ++i) {
        coeff[i] = new DataType **[4]; // [a, b, aa, bb]
        for (int j = 0; j < 4; ++j) {
            coeff[i][j] = new DataType *[3]; // [i-1, i, i+1]
            for (int k = 0; k < 3; ++k) {
                coeff[i][j][k] = new DataType[n];
            }
        }
    }

    DataType a1, b1, b2, b3;
    a1 = 0.875;
    b1 = 0.1251282341599089;
    b2 = -0.2487176584009104;
    b3 = 0.0001282341599089;

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {

        coeff[0][0][0][i] = a1;
        coeff[0][0][1][i] = 1.0;
        coeff[0][0][2][i] = 0.0;
        coeff[1][0][0][i] = 0.0;
        coeff[1][0][1][i] = 1.0;
        coeff[1][0][2][i] = a1;

        coeff[0][1][0][i] = b1 * h;
        coeff[0][1][1][i] = b2 * h;
        coeff[0][1][2][i] = b3 * h;
        coeff[1][1][0][i] = -b3 * h;
        coeff[1][1][1][i] = -b2 * h;
        coeff[1][1][2][i] = -b1 * h;

        coeff[0][2][0][i] = -9.0 / 8.0 / h;
        coeff[0][2][1][i] = 0.0;
        coeff[0][2][2][i] = 9.0 / 8.0 / h;
        coeff[1][2][0][i] = -9.0 / 8.0 / h;
        coeff[1][2][1][i] = 0.0;
        coeff[1][2][2][i] = 9.0 / 8.0 / h;

        coeff[0][3][0][i] = -1.0 / 8.0;
        coeff[0][3][1][i] = 1.0;
        coeff[0][3][2][i] = -1.0 / 8.0;
        coeff[1][3][0][i] = -1.0 / 8.0;
        coeff[1][3][1][i] = 1.0;
        coeff[1][3][2][i] = -1.0 / 8.0;

    }

    ccd_coefficient_boundary_condition(coeff[0], n, h);
    ccd_coefficient_boundary_condition(coeff[1], n, h);

    twin_dec(coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], n);
    twin_dec(coeff[1][0], coeff[1][1], coeff[1][2], coeff[1][3], n);

    return coeff;
}


void uccd_find_fx(scalar_data *f, structured_grid *geo, vector_data *vel) {

    DataType c1, c2, c3;
    c3 = -0.06096119008109;
    c2 = 1.99692238016218;
    c1 = -1.93596119008109;

#pragma omp parallel for
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto coeff = uccd_coefficient_matrix(f->Nx, geo->dx);

            DataType **su, **sd;
            su = new DataType *[2];
            sd = new DataType *[2];
            for (int i = 0; i < 2; ++i) {
                su[i] = new DataType[f->Nx];
                sd[i] = new DataType[f->Nx];
            }

            su[0][0] = (-3.5 * f->flux[0][j][k]
                        + 4.0 * f->flux[1][j][k]
                        - 0.5 * f->flux[2][j][k]) / geo->dx;
            su[1][0] = (34.0 / 3.0 * f->flux[0][j][k]
                        - 83.0 / 4.0 * f->flux[1][j][k]
                        + 10.0 * f->flux[2][j][k]
                        - 7.0 / 12.0 * f->flux[3][j][k]) / geo->dx / geo->dx;

            su[0][f->Nx - 1] = -(-3.5 * f->flux[f->Nx - 1][j][k]
                                 + 4.0 * f->flux[f->Nx - 2][j][k]
                                 - 0.5 * f->flux[f->Nx - 3][j][k]) / geo->dx;
            su[1][f->Nx - 1] = (34.0 / 3.0 * f->flux[f->Nx - 1][j][k]
                                - 83.0 / 4.0 * f->flux[f->Nx - 2][j][k]
                                + 10.0 * f->flux[f->Nx - 3][j][k]
                                - 7.0 / 12.0 * f->flux[f->Nx - 4][j][k]) / geo->dx / geo->dx;

            sd[0][0] = su[0][0];
            sd[1][0] = su[1][0];
            sd[0][f->Nx - 1] = su[0][f->Nx - 1];
            sd[1][f->Nx - 1] = su[1][f->Nx - 1];

            for (int i = 1; i < f->Nx - 1; ++i) {
                su[0][i] = (c1 * f->flux[i - 1][j][k]
                            + c2 * f->flux[i][j][k]
                            + c3 * f->flux[i + 1][j][k]) / geo->dx;
                su[1][i] = (3.0 * f->flux[i - 1][j][k] - 6.0 * f->flux[i][j][k]
                            + 3.0 * f->flux[i + 1][j][k]) / geo->dx / geo->dx;

                sd[0][i] = -(c3 * f->flux[i - 1][j][k]
                             + c2 * f->flux[i][j][k]
                             + c1 * f->flux[i + 1][j][k]) / geo->dx;
                sd[1][i] = su[1][i];
            }

            twin_bks(coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], su[0], su[1], f->Nx);
            twin_bks(coeff[1][0], coeff[1][1], coeff[1][2], coeff[1][3], sd[0], sd[1], f->Nx);

            for (int i = 0; i < f->Nx; ++i) {
                if (vel->x.data[i][j][k] > 0.0) {
                    f->fx[i][j][k] = su[0][i];
                } else {
                    f->fx[i][j][k] = sd[0][i];
                }
                f->fxx[i][j][k] = 0.5 * (su[1][i] + sd[1][i]);
            }

            delete4d(coeff, 2, 4, 3);
            delete2d(su, 2);
            delete2d(sd, 2);

        }
    }
}

void uccd_find_fy(scalar_data *f, structured_grid *geo, vector_data *vel) {

    DataType c1, c2, c3;
    c3 = -0.06096119008109;
    c2 = 1.99692238016218;
    c1 = -1.93596119008109;

#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int k = 0; k < f->Nz; ++k) {
            auto coeff = uccd_coefficient_matrix(f->Ny, geo->dy);

            DataType **su, **sd;
            su = new DataType *[2];
            sd = new DataType *[2];
            for (int i = 0; i < 2; ++i) {
                su[i] = new DataType[f->Ny];
                sd[i] = new DataType[f->Ny];
            }

            su[0][0] = (-3.5 * f->flux[i][0][k]
                        + 4.0 * f->flux[i][1][k]
                        - 0.5 * f->flux[i][2][k]) / geo->dy;
            su[1][0] = (34.0 / 3.0 * f->flux[i][0][k]
                        - 83.0 / 4.0 * f->flux[i][1][k]
                        + 10.0 * f->flux[i][2][k]
                        - 7.0 / 12.0 * f->flux[i][3][k]) / geo->dy / geo->dy;

            su[0][f->Ny - 1] = -(-3.5 * f->flux[i][f->Ny - 1][k]
                                 + 4.0 * f->flux[i][f->Ny - 2][k]
                                 - 0.5 * f->flux[i][f->Ny - 3][k]) / geo->dy;
            su[1][f->Ny - 1] = (34.0 / 3.0 * f->flux[i][f->Ny - 1][k]
                                - 83.0 / 4.0 * f->flux[i][f->Ny - 2][k]
                                + 10.0 * f->flux[i][f->Ny - 3][k]
                                - 7.0 / 12.0 * f->flux[i][f->Ny - 4][k]) / geo->dy / geo->dy;

            sd[0][0] = su[0][0];
            sd[1][0] = su[1][0];

            sd[0][f->Ny - 1] = su[0][f->Ny - 1];
            sd[1][f->Ny - 1] = su[1][f->Ny - 1];

            for (int j = 1; j < f->Ny - 1; ++j) {
                su[0][j] = (c1 * f->flux[i][j - 1][k]
                            + c2 * f->flux[i][j][k]
                            + c3 * f->flux[i][j + 1][k]) / geo->dy;
                su[1][j] = (3.0 * f->flux[i][j - 1][k] - 6.0 * f->flux[i][j][k]
                            + 3.0 * f->flux[i][j + 1][k]) / geo->dy / geo->dy;

                sd[0][j] = -(c3 * f->flux[i][j - 1][k]
                             + c2 * f->flux[i][j][k]
                             + c1 * f->flux[i][j + 1][k]) / geo->dy;
                sd[1][j] = su[1][j];
            }

            twin_bks(coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], su[0], su[1], f->Ny);
            twin_bks(coeff[1][0], coeff[1][1], coeff[1][2], coeff[1][3], sd[0], sd[1], f->Ny);

            for (int j = 0; j < f->Ny; ++j) {
                if (vel->y.data[i][j][k] > 0.0) {
                    f->fy[i][j][k] = su[0][j];
                } else {
                    f->fy[i][j][k] = sd[0][j];
                }
                f->fyy[i][j][k] = 0.5 * (su[1][j] + sd[1][j]);
            }

            delete4d(coeff, 2, 4, 3);
            delete2d(su, 2);
            delete2d(sd, 2);

        }
    }
}

void uccd_find_fz(scalar_data *f, structured_grid *geo, vector_data *vel) {

    DataType c1, c2, c3;
    c3 = -0.06096119008109;
    c2 = 1.99692238016218;
    c1 = -1.93596119008109;

#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            auto coeff = uccd_coefficient_matrix(f->Nz, geo->dz);

            DataType **su, **sd;
            su = new DataType *[2];
            sd = new DataType *[2];
            for (int i = 0; i < 2; ++i) {
                su[i] = new DataType[f->Nz];
                sd[i] = new DataType[f->Nz];
            }

            su[0][0] = (-3.5 * f->flux[i][j][0]
                        + 4.0 * f->flux[i][j][1]
                        - 0.5 * f->flux[i][j][2]) / geo->dz;
            su[1][0] = (34.0 / 3.0 * f->flux[i][j][0]
                        - 83.0 / 4.0 * f->flux[i][j][1]
                        + 10.0 * f->flux[i][j][2]
                        - 7.0 / 12.0 * f->flux[i][j][3]) / geo->dz / geo->dz;

            su[0][f->Nz - 1] = -(-3.5 * f->flux[i][j][f->Nz - 1]
                                 + 4.0 * f->flux[i][j][f->Nz - 2]
                                 - 0.5 * f->flux[i][j][f->Nz - 3]) / geo->dz;
            su[1][f->Nz - 1] = (34.0 / 3.0 * f->flux[i][j][f->Nz - 1]
                                - 83.0 / 4.0 * f->flux[i][j][f->Nz - 2]
                                + 10.0 * f->flux[i][j][f->Nz - 3]
                                - 7.0 / 12.0 * f->flux[i][j][f->Nz - 4]) / geo->dz / geo->dz;

            sd[0][0] = su[0][0];
            sd[1][0] = su[1][0];

            sd[0][f->Nz - 1] = su[0][f->Nz - 1];
            sd[1][f->Nz - 1] = su[1][f->Nz - 1];

            for (int k = 1; k < f->Nz - 1; ++k) {
                su[0][k] = (c1 * f->flux[i][j][k - 1]
                            + c2 * f->flux[i][j][k]
                            + c3 * f->flux[i][j][k + 1]) / geo->dz;
                su[1][k] = (3.0 * f->flux[i][j][k - 1] - 6.0 * f->flux[i][j][k]
                            + 3.0 * f->flux[i][j][k + 1]) / geo->dz / geo->dz;

                sd[0][k] = -(c3 * f->flux[i][j][k - 1]
                             + c2 * f->flux[i][j][k]
                             + c1 * f->flux[i][j][k + 1]) / geo->dz;
                sd[1][k] = su[1][k];
            }

            twin_bks(coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], su[0], su[1], f->Nz);
            twin_bks(coeff[1][0], coeff[1][1], coeff[1][2], coeff[1][3], sd[0], sd[1], f->Nz);

            for (int k = 0; k < f->Nz; ++k) {
                if (vel->z.data[i][j][k] > 0.0) {
                    f->fz[i][j][k] = su[0][k];
                } else {
                    f->fz[i][j][k] = sd[0][k];
                }
                f->fzz[i][j][k] = 0.5 * (su[1][k] + sd[1][k]);
            }

            delete4d(coeff, 2, 4, 3);
            delete2d(su, 2);
            delete2d(sd, 2);

        }
    }

}

void
uccd_find_derivatives(scalar_data *f, structured_grid *geo, vector_data *vel) {
    uccd_find_fx(f, geo, vel);
    if (f->ndim > 1) {
        uccd_find_fy(f, geo, vel);
    }
    if (f->ndim > 2) {
        uccd_find_fz(f, geo, vel);
    }
}

