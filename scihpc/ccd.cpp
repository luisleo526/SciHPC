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
    for (int j = 1; j <= f->ny; ++j) {
        for (int k = 1; k <= f->nz; ++k) {

            auto index = f->index_mapping(1, j, k);
            auto coeff = ccd_coefficient_matrix(f->Nx, geo->dx);

            // Find src
            DataType *s, *ss;
            s = new DataType[f->Nx];
            ss = new DataType[f->Nx];

            // Src boundary condition
            s[0] = (-3.5 * f->flux[index.i][index.j][index.k]
                    + 4.0 * f->flux[index.i + 1][index.j][index.k]
                    - 0.5 * f->flux[index.i + 2][index.j][index.k]) / geo->dx;

            ss[0] = (34.0 / 3.0 * f->flux[index.i][index.j][index.k]
                     - 83.0 / 4.0 * f->flux[index.i + 1][index.j][index.k]
                     + 10.0 * f->flux[index.i + 2][index.j][index.k]
                     - 7.0 / 12.0 * f->flux[index.i + 3][index.j][index.k]) / geo->dx / geo->dx;

            s[f->Nx - 1] = -(-3.5 * f->flux[index.i + f->Nx - 1][index.j][index.k]
                             + 4.0 * f->flux[index.i + f->Nx - 2][index.j][index.k]
                             - 0.5 * f->flux[index.i + f->Nx - 3][index.j][index.k]) / geo->dx;

            ss[f->Nx - 1] = (34.0 / 3.0 * f->flux[index.i + f->Nx - 1][index.j][index.k]
                             - 83.0 / 4.0 * f->flux[index.i + f->Nx - 2][index.j][index.k]
                             + 10.0 * f->flux[index.i + f->Nx - 3][index.j][index.k]
                             - 7.0 / 12.0 * f->flux[index.i + f->Nx - 4][index.j][index.k]) / geo->dx / geo->dx;

            for (int i = 1; i < f->Nx - 1; ++i) {
                s[i] = 15.0 / 16.0 * (f->flux[index.i + i + 1][index.j][index.k]
                                      - f->flux[index.i + i - 1][index.j][index.k]) / geo->dx;
                ss[i] = (3.0 * f->flux[index.i + i - 1][index.j][index.k]
                         - 6.0 * f->flux[index.i + i][index.j][index.k]
                         + 3.0 * f->flux[index.i + i + 1][index.j][index.k]) / geo->dx / geo->dx;
            }

            twin_bks(coeff[0], coeff[1], coeff[2], coeff[3], s, ss, f->Nx);
            for (int i = 0; i < f->Nx; ++i) {
                f->fx[index.i + i][index.j][index.k] = s[i];
                f->fxx[index.i + i][index.j][index.k] = ss[i];
            }

            delete3d(coeff, 4, 3);
            delete[] s;
            delete[] ss;
        }
    }
}

void ccd_find_fy(scalar_data *f, structured_grid *geo) {

#pragma omp parallel for
    for (int i = 1; i <= f->nx; ++i) {
        for (int k = 1; k <= f->nz; ++k) {

            auto index = f->index_mapping(i, 1, k);
            auto coeff = ccd_coefficient_matrix(f->Ny, geo->dy);

            // Find src
            DataType *s, *ss;
            s = new DataType[f->Ny];
            ss = new DataType[f->Ny];

            s[0] = (-3.5 * f->flux[index.i][index.j][index.k]
                    + 4.0 * f->flux[index.i][index.j + 1][index.k]
                    - 0.5 * f->flux[index.i][index.j + 2][index.k]) / geo->dy;

            ss[0] = (34.0 / 3.0 * f->flux[index.i][index.j][index.k]
                     - 83.0 / 4.0 * f->flux[index.i][index.j + 1][index.k]
                     + 10.0 * f->flux[index.i][index.j + 2][index.k]
                     - 7.0 / 12.0 * f->flux[index.i][index.j + 3][index.k]) / geo->dy / geo->dy;

            s[f->Ny - 1] = -(-3.5 * f->flux[index.i][index.j + f->Ny - 1][index.k]
                             + 4.0 * f->flux[index.i][index.j + f->Ny - 2][index.k]
                             - 0.5 * f->flux[index.i][index.j + f->Ny - 3][index.k]) / geo->dy;

            ss[f->Ny - 1] = (34.0 / 3.0 * f->flux[index.i][index.j + f->Ny - 1][index.k]
                             - 83.0 / 4.0 * f->flux[index.i][index.j + f->Ny - 2][index.k]
                             + 10.0 * f->flux[index.i][index.j + f->Ny - 3][index.k]
                             - 7.0 / 12.0 * f->flux[index.i][index.j + f->Ny - 4][index.k]) / geo->dy / geo->dy;

            for (int j = 1; j < f->Ny - 1; ++j) {
                s[j] = 15.0 / 16.0 * (f->flux[index.i][index.j + j + 1][index.k]
                                      - f->flux[index.i][index.j + j - 1][index.k]) / geo->dy;
                ss[j] = (3.0 * f->flux[index.i][index.j + j - 1][index.k]
                         - 6.0 * f->flux[index.i][index.j + j][index.k]
                         + 3.0 * f->flux[index.i][index.j + j + 1][index.k]) / geo->dy / geo->dy;
            }

            twin_bks(coeff[0], coeff[1], coeff[2], coeff[3], s, ss, f->Ny);
            for (int j = 0; j < f->Ny; ++j) {
                f->fy[index.i][index.j + j][index.k] = s[j];
                f->fyy[index.i][index.j + j][index.k] = ss[j];
            }

            delete3d(coeff, 4, 3);
            delete[] s;
            delete[] ss;
        }
    }
}

void ccd_find_fz(scalar_data *f, structured_grid *geo) {

#pragma omp parallel for
    for (int i = 1; i <= f->nx; ++i) {
        for (int j = 1; j <= f->ny; ++j) {

            auto index = f->index_mapping(i, j, 1);
            auto coeff = ccd_coefficient_matrix(f->Nz, geo->dz);

            // Find src
            DataType *s, *ss;
            s = new DataType[f->Nz];
            ss = new DataType[f->Nz];

            s[0] = (-3.5 * f->flux[index.i][index.j][index.k]
                    + 4.0 * f->flux[index.i][index.j][index.k + 1]
                    - 0.5 * f->flux[index.i][index.j][index.k + 2]) / geo->dz;

            ss[0] = (34.0 / 3.0 * f->flux[index.i][index.j][index.k]
                     - 83.0 / 4.0 * f->flux[index.i][index.j][index.k + 1]
                     + 10.0 * f->flux[index.i][index.j][index.k + 2]
                     - 7.0 / 12.0 * f->flux[index.i][index.j][index.k + 3]) / geo->dz / geo->dz;

            s[f->Nz - 1] = -(-3.5 * f->flux[index.i][index.j][index.k + f->Nz - 1]
                             + 4.0 * f->flux[index.i][index.j][index.k + f->Nz - 2]
                             - 0.5 * f->flux[index.i][index.j][index.k + f->Nz - 3]) / geo->dz;

            ss[f->Nz - 1] = (34.0 / 3.0 * f->flux[index.i][index.j][index.k + f->Nz - 1]
                             - 83.0 / 4.0 * f->flux[index.i][index.j][index.k + f->Nz - 2]
                             + 10.0 * f->flux[index.i][index.j][index.k + f->Nz - 3]
                             - 7.0 / 12.0 * f->flux[index.i][index.j][index.k + f->Nz - 4]) / geo->dz / geo->dz;

            for (int k = 1; k < f->Nz - 1; ++k) {
                s[k] = 15.0 / 16.0 * (f->flux[index.i][index.j][index.k + k + 1]
                                      - f->flux[index.i][index.j][index.k + k - 1]) / geo->dz;
                ss[k] = (3.0 * f->flux[index.i][index.j][index.k + k - 1]
                         - 6.0 * f->flux[index.i][index.j][index.k + k]
                         + 3.0 * f->flux[index.i][index.j][index.k + k + 1]) / geo->dz / geo->dz;
            }

            twin_bks(coeff[0], coeff[1], coeff[2], coeff[3], s, ss, f->Nz);

            for (int k = 0; k < f->Nz; ++k) {
                f->fz[index.i][index.j][index.k + k] = s[k];
                f->fzz[index.i][index.j][index.k + k] = ss[k];
            }

            delete3d(coeff, 4, 3);
            delete[] s;
            delete[] ss;

        }
    }
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
    for (int j = 1; j <= f->ny; ++j) {
        for (int k = 1; k <= f->nz; ++k) {
            auto index = f->index_mapping(-2, j, k);
            auto coeff = uccd_coefficient_matrix(f->Nx, geo->dx);

            DataType **su, **sd;
            su = new DataType *[2];
            sd = new DataType *[2];
            for (int i = 0; i < 2; ++i) {
                su[i] = new DataType[f->Nx];
                sd[i] = new DataType[f->Nx];
            }

            su[0][0] = (-3.5 * f->flux[index.i][index.j][index.k]
                        + 4.0 * f->flux[index.i + 1][index.j][index.k]
                        - 0.5 * f->flux[index.i + 2][index.j][index.k]) / geo->dx;
            sd[0][0] = su[0][0];

            su[1][0] = (34.0 / 3.0 * f->flux[index.i][index.j][index.k]
                        - 83.0 / 4.0 * f->flux[index.i + 1][index.j][index.k]
                        + 10.0 * f->flux[index.i + 2][index.j][index.k]
                        - 7.0 / 12.0 * f->flux[index.i + 3][index.j][index.k]) / geo->dx / geo->dx;
            sd[1][0] = su[1][0];

            su[0][f->Nx - 1] = -(-3.5 * f->flux[index.i + f->Nx - 1][index.j][index.k]
                                 + 4.0 * f->flux[index.i + f->Nx - 2][index.j][index.k]
                                 - 0.5 * f->flux[index.i + f->Nx - 3][index.j][index.k]) / geo->dx;
            sd[0][f->Nx - 1] = su[0][f->Nx - 1];

            su[1][f->Nx - 1] = (34.0 / 3.0 * f->flux[index.i + f->Nx - 1][index.j][index.k]
                                - 83.0 / 4.0 * f->flux[index.i + f->Nx - 2][index.j][index.k]
                                + 10.0 * f->flux[index.i + f->Nx - 3][index.j][index.k]
                                - 7.0 / 12.0 * f->flux[index.i + f->Nx - 4][index.j][index.k]) / geo->dx / geo->dx;
            sd[1][f->Nx - 1] = su[1][f->Nx - 1];

            for (int i = 1; i < f->Nx - 1; ++i) {
                su[0][i] = (c1 * f->flux[index.i + i - 1][index.j][index.k]
                            + c2 * f->flux[index.i + i][index.j][index.k]
                            + c3 * f->flux[index.i + i + 1][index.j][index.k]) / geo->dx;

                sd[0][i] = -(c3 * f->flux[index.i + i - 1][index.j][index.k]
                             + c2 * f->flux[index.i + i][index.j][index.k]
                             + c1 * f->flux[index.i + i + 1][index.j][index.k]) / geo->dx;

                su[1][i] = (3.0 * f->flux[index.i + i - 1][index.j][index.k]
                            - 6.0 * f->flux[index.i + i][index.j][index.k]
                            + 3.0 * f->flux[index.i + i + 1][index.j][index.k]) / geo->dx / geo->dx;
                sd[1][i] = su[1][i];
            }

            twin_bks(coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], su[0], su[1], f->Nx);
            twin_bks(coeff[1][0], coeff[1][1], coeff[1][2], coeff[1][3], sd[0], sd[1], f->Nx);

            for (int i = 0; i < f->Nx; ++i) {
                if (vel->x.data[index.i + i][index.j][index.k] > 0.0) {
                    f->fx[index.i + i][index.j][index.k] = su[0][i];
                } else {
                    f->fx[index.i + i][index.j][index.k] = sd[0][i];
                }
                f->fxx[index.i + i][index.j][index.k] = 0.5 * (su[1][i] + sd[1][i]);
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
    for (int i = 1; i <= f->nx; ++i) {
        for (int k = 1; k <= f->nz; ++k) {
            auto index = f->index_mapping(i, 1, k);
            auto coeff = uccd_coefficient_matrix(f->Ny, geo->dy);

            DataType **su, **sd;
            su = new DataType *[2];
            sd = new DataType *[2];
            for (int j = 0; j < 2; ++j) {
                su[j] = new DataType[f->Ny];
                sd[j] = new DataType[f->Ny];
            }

            su[0][0] = (-3.5 * f->flux[index.i][index.j][index.k]
                        + 4.0 * f->flux[index.i][index.j + 1][index.k]
                        - 0.5 * f->flux[index.i][index.j + 2][index.k]) / geo->dy;
            sd[0][0] = su[0][0];

            su[1][0] = (34.0 / 3.0 * f->flux[index.i][index.j][index.k]
                        - 83.0 / 4.0 * f->flux[index.i][index.j + 1][index.k]
                        + 10.0 * f->flux[index.i][index.j + 2][index.k]
                        - 7.0 / 12.0 * f->flux[index.i][index.j + 3][index.k]) / geo->dy / geo->dy;
            sd[1][0] = su[1][0];

            su[0][f->Ny - 1] = -(-3.5 * f->flux[index.i][index.j + f->Ny - 1][index.k]
                                 + 4.0 * f->flux[index.i][index.j + f->Ny - 2][index.k]
                                 - 0.5 * f->flux[index.i][index.j + f->Ny - 3][index.k]) / geo->dy;
            sd[0][f->Ny - 1] = su[0][f->Ny - 1];

            su[1][f->Ny - 1] = (34.0 / 3.0 * f->flux[index.i][index.j + f->Ny - 1][index.k]
                                - 83.0 / 4.0 * f->flux[index.i][index.j + f->Ny - 2][index.k]
                                + 10.0 * f->flux[index.i][index.j + f->Ny - 3][index.k]
                                - 7.0 / 12.0 * f->flux[index.i][index.j + f->Ny - 4][index.k]) / geo->dy / geo->dy;
            sd[1][f->Ny - 1] = su[1][f->Ny - 1];

            for (int j = 1; j < f->Ny - 1; ++j) {
                su[0][j] = (c1 * f->flux[index.i][index.j + j - 1][index.k]
                            + c2 * f->flux[index.i][index.j + j][index.k]
                            + c3 * f->flux[index.i][index.j + j + 1][index.k]) / geo->dy;

                sd[0][j] = -(c3 * f->flux[index.i][index.j + j - 1][index.k]
                             + c2 * f->flux[index.i][index.j + j][index.k]
                             + c1 * f->flux[index.i][index.j + j + 1][index.k]) / geo->dy;

                su[1][j] = (3.0 * f->flux[index.i][index.j + j - 1][index.k]
                            - 6.0 * f->flux[index.i][index.j + j][index.k]
                            + 3.0 * f->flux[index.i][index.j + j + 1][index.k]) / geo->dy / geo->dy;
                sd[1][j] = su[1][j];
            }

            twin_bks(coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], su[0], su[1], f->Ny);
            twin_bks(coeff[1][0], coeff[1][1], coeff[1][2], coeff[1][3], sd[0], sd[1], f->Ny);

            for (int j = 0; j < f->Ny; ++j) {
                if (vel->y.data[index.i][index.j + j][index.k] > 0.0) {
                    f->fy[index.i][index.j + j][index.k] = su[0][j];
                } else {
                    f->fy[index.i][index.j + j][index.k] = sd[0][j];
                }
                f->fyy[index.i][index.j + j][index.k] = 0.5 * (su[1][j] + sd[1][j]);
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
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            auto index = f->index_mapping(i, j, 1);
            auto coeff = uccd_coefficient_matrix(f->Nz, geo->dz);

            DataType **su, **sd;
            su = new DataType *[2];
            sd = new DataType *[2];
            for (int k = 0; k < 2; ++k) {
                su[k] = new DataType[f->Nz];
                sd[k] = new DataType[f->Nz];
            }

            su[0][0] = (-3.5 * f->flux[index.i][index.j][index.k]
                        + 4.0 * f->flux[index.i][index.j][index.k + 1]
                        - 0.5 * f->flux[index.i][index.j][index.k + 2]) / geo->dz;
            sd[0][0] = su[0][0];

            su[1][0] = (34.0 / 3.0 * f->flux[index.i][index.j][index.k]
                        - 83.0 / 4.0 * f->flux[index.i][index.j][index.k + 1]
                        + 10.0 * f->flux[index.i][index.j][index.k + 2]
                        - 7.0 / 12.0 * f->flux[index.i][index.j][index.k + 3]) / geo->dz / geo->dz;
            sd[1][0] = su[1][0];

            su[0][f->Nz - 1] = -(-3.5 * f->flux[index.i][index.j][index.k + f->Nz - 1]
                                 + 4.0 * f->flux[index.i][index.j][index.k + f->Nz - 2]
                                 - 0.5 * f->flux[index.i][index.j][index.k + f->Nz - 3]) / geo->dz;
            sd[0][f->Nz - 1] = su[0][f->Nz - 1];

            su[1][f->Nz - 1] = (34.0 / 3.0 * f->flux[index.i][index.j][index.k + f->Nz - 1]
                                - 83.0 / 4.0 * f->flux[index.i][index.j][index.k + f->Nz - 2]
                                + 10.0 * f->flux[index.i][index.j][index.k + f->Nz - 3]
                                - 7.0 / 12.0 * f->flux[index.i][index.j][index.k + f->Nz - 4]) / geo->dz / geo->dz;

            sd[1][f->Nz - 1] = su[1][f->Nz - 1];

            for (int k = 1; k < f->Nz - 1; ++k) {
                su[0][k] = (c1 * f->flux[index.i][index.j][index.k + k - 1]
                            + c2 * f->flux[index.i][index.j][index.k + k]
                            + c3 * f->flux[index.i][index.j][index.k + k + 1]) / geo->dz;

                sd[0][k] = -(c3 * f->flux[index.i][index.j][index.k + k - 1]
                             + c2 * f->flux[index.i][index.j][index.k + k]
                             + c1 * f->flux[index.i][index.j][index.k + k + 1]) / geo->dz;

                su[1][k] = (3.0 * f->flux[index.i][index.j][index.k + k - 1]
                            - 6.0 * f->flux[index.i][index.j][index.k + k]
                            + 3.0 * f->flux[index.i][index.j][index.k + k + 1]) / geo->dz / geo->dz;
                sd[1][k] = su[1][k];
            }

            twin_bks(coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], su[0], su[1], f->Nz);
            twin_bks(coeff[1][0], coeff[1][1], coeff[1][2], coeff[1][3], sd[0], sd[1], f->Nz);

            for (int k = 0; k < f->Nz; ++k) {
                if (vel->z.data[index.i][index.j][index.k + k] > 0.0) {
                    f->fz[index.i][index.j][index.k + k] = su[0][k];
                } else {
                    f->fz[index.i][index.j][index.k + k] = sd[0][k];
                }
                f->fzz[index.i][index.j][index.k + k] = 0.5 * (su[1][k] + sd[1][k]);
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

