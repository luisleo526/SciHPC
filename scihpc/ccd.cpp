//
// Created by 溫晧良 on 2022/9/30.
//

#include "ccd.h"

void ccd_coefficient_boundary_condition(DataType ***coeff, const int n, const DataType dx) {

    coeff[0][1][0] = 1.0;
    coeff[1][1][0] = 0.0;
    coeff[2][1][0] = 0.0;
    coeff[3][1][0] = 1.0;

    coeff[0][2][0] = 2.0;
    coeff[1][2][0] = -dx;
    coeff[2][2][0] = -2.5 / dx;
    coeff[3][2][0] = 8.5;

    coeff[0][1][n - 1] = 1.0;
    coeff[1][1][n - 1] = 0.0;
    coeff[2][1][n - 1] = 0.0;
    coeff[3][1][n - 1] = 1.0;

    coeff[0][0][n - 1] = 2.0;
    coeff[1][0][n - 1] = dx;
    coeff[2][0][n - 1] = 2.5 / dx;
    coeff[3][0][n - 1] = 8.5;

}

DataType ***ccd_coefficient_matrix(const int n, const DataType dx) {

    auto coeff = new DataType **[4];
    for (int i = 0; i < 4; ++i) {
        coeff[i] = new DataType *[3];
        for (int j = 0; j < 3; ++j) {
            coeff[i][j] = new DataType[n];
        }
    }

    for (int i = 0; i < n; ++i) {
        coeff[0][0][i] = 7.0 / 16.0;
        coeff[0][1][i] = 1.0;
        coeff[0][2][i] = 7.0 / 16.0;

        coeff[1][0][i] = dx / 16.0;
        coeff[1][1][i] = 0.0;
        coeff[1][2][i] = -dx / 16.0;

        coeff[2][0][i] = -9.0 / 8.0 / dx;
        coeff[2][1][i] = 0.0;
        coeff[2][2][i] = 9.0 / 8.0 / dx;

        coeff[3][0][i] = -1.0 / 8.0;
        coeff[3][1][i] = 1.0;
        coeff[3][2][i] = -1.0 / 8.0;
    }

    ccd_coefficient_boundary_condition(coeff, n, dx);
    twin_dec(coeff[0], coeff[1], coeff[2], coeff[3], n);

    return coeff;
}

void ccd_find_fx(scalar_data *f, const DataType dx) {
    for (int j = 1; j <= f->ny; ++j) {
        for (int k = 1; k <= f->nz; ++k) {

            auto index = f->index_mapping(1, j, k);
            auto coeff = ccd_coefficient_matrix(f->nx, dx);

            // Find src
            DataType *s, *ss;
            s = new DataType[f->nx];
            ss = new DataType[f->nx];

            s[0] = (-3.5 * f->data[index.i][index.j][index.k]
                    + 4.0 * f->data[index.i + 1][index.j][index.k]
                    - 0.5 * f->data[index.i + 2][index.j][index.k]) / dx;

            ss[0] = (34.0 / 3.0 * f->data[index.i][index.j][index.k]
                     - 83.0 / 4.0 * f->data[index.i + 1][index.j][index.k]
                     + 10.0 * f->data[index.i + 2][index.j][index.k]
                     - 7.0 / 12.0 * f->data[index.i + 3][index.j][index.k]) / dx / dx;

            s[f->nx - 1] = -(-3.5 * f->data[index.i + f->nx - 1][index.j][index.k]
                             + 4.0 * f->data[index.i + f->nx - 2][index.j][index.k]
                             - 0.5 * f->data[index.i + f->nx - 3][index.j][index.k]) / dx;

            ss[f->nx - 1] = (34.0 / 3.0 * f->data[index.i + f->nx - 1][index.j][index.k]
                             - 83.0 / 4.0 * f->data[index.i + f->nx - 2][index.j][index.k]
                             + 10.0 * f->data[index.i + f->nx - 3][index.j][index.k]
                             - 7.0 / 12.0 * f->data[index.i + f->nx - 4][index.j][index.k]) / dx / dx;

            for (int i = 1; i < f->nx - 1; ++i) {
                s[i] = 15.0 / 16.0 * (f->data[index.i + i + 1][index.j][index.k]
                                      - f->data[index.i + i - 1][index.j][index.k]) / dx;
                ss[i] = (3.0 * f->data[index.i + i - 1][index.j][index.k]
                         - 6.0 * f->data[index.i + i][index.j][index.k]
                         + 3.0 * f->data[index.i + i + 1][index.j][index.k]) / dx / dx;
            }

            twin_bks(coeff[0], coeff[1], coeff[2], coeff[3], s, ss, f->nx);
            for (int i = 0; i < f->nx; ++i) {
                index = f->index_mapping(i + 1, j, k);
                f->fx[index.i][index.j][index.k] = s[i];
                f->fxx[index.i][index.j][index.k] = ss[i];
            }

            delete3d(coeff, 4, 3);
            delete[] s;
            delete[] ss;
        }
    }
}

void ccd_find_fy(scalar_data *f, DataType dy) {
    for (int i = 1; i <= f->nx; ++i) {
        for (int k = 1; k <= f->nz; ++k) {

            auto index = f->index_mapping(i, 1, k);
            auto coeff = ccd_coefficient_matrix(f->ny, dy);

            // Find src
            DataType *s, *ss;
            s = new DataType[f->ny];
            ss = new DataType[f->ny];

            s[0] = (-3.5 * f->data[index.i][index.j][index.k]
                    + 4.0 * f->data[index.i][index.j + 1][index.k]
                    - 0.5 * f->data[index.i][index.j + 2][index.k]) / dy;

            ss[0] = (34.0 / 3.0 * f->data[index.i][index.j][index.k]
                     - 83.0 / 4.0 * f->data[index.i][index.j + 1][index.k]
                     + 10.0 * f->data[index.i][index.j + 2][index.k]
                     - 7.0 / 12.0 * f->data[index.i][index.j + 3][index.k]) / dy / dy;

            s[f->ny - 1] = -(-3.5 * f->data[index.i][index.j + f->ny - 1][index.k]
                             + 4.0 * f->data[index.i][index.j + f->ny - 2][index.k]
                             - 0.5 * f->data[index.i][index.j + f->ny - 3][index.k]) / dy;

            ss[f->ny - 1] = (34.0 / 3.0 * f->data[index.i][index.j + f->ny - 1][index.k]
                             - 83.0 / 4.0 * f->data[index.i][index.j + f->ny - 2][index.k]
                             + 10.0 * f->data[index.i][index.j + f->ny - 3][index.k]
                             - 7.0 / 12.0 * f->data[index.i][index.j + f->ny - 4][index.k]) / dy / dy;

            for (int j = 1; j < f->ny - 1; ++j) {
                s[j] = 15.0 / 16.0 * (f->data[index.i][index.j + j + 1][index.k]
                                      - f->data[index.i][index.j + j - 1][index.k]) / dy;
                ss[j] = (3.0 * f->data[index.i][index.j + j - 1][index.k]
                         - 6.0 * f->data[index.i][index.j + j][index.k]
                         + 3.0 * f->data[index.i][index.j + j + 1][index.k]) / dy / dy;
            }

            twin_bks(coeff[0], coeff[1], coeff[2], coeff[3], s, ss, f->ny);
            for (int j = 0; j < f->ny; ++j) {
                index = f->index_mapping(i, j + 1, k);
                f->fy[index.i][index.j][index.k] = s[j];
                f->fyy[index.i][index.j][index.k] = ss[j];
            }

            delete3d(coeff, 4, 3);
            delete[] s;
            delete[] ss;
        }
    }
}

