//
// Created by leo on 10/3/22.
//

#include "weno_solver.h"

weno_solver::weno_solver(scalar_data *f, structured_grid *geo) {
    fp = init_array(f->Nx, f->Ny, f->Nz);
    fm = init_array(f->Nx, f->Ny, f->Nz);
    fh = init_array(f->Nx, f->Ny, f->Nz);
    dx = geo->dx;
    dy = geo->dy;
    dz = geo->dz;
}

weno_data weno_solver::indicators_p(DataType a, DataType b, DataType c, DataType d, DataType e) {

    auto b0 = 13.0 * (c - 2.0 * d + e) * (c - 2.0 * d + e) + 3.0 * (3.0 * c - 4.0 * d + e) * (3.0 * c - 4.0 * d + e);
    auto b1 = 13.0 * (b - 2.0 * c + d) * (b - 2.0 * c + d) + 3.0 * (b - d) * (b - d);
    auto b2 = 13.0 * (a - 2.0 * b + c) * (a - 2.0 * b + c) + 3.0 * (a - 4.0 * b + 3.0 * c) * (a - 4.0 * b + 3.0 * c);

    return weno_data{b0, b1, b2};
}

weno_data weno_solver::indicators_m(DataType a, DataType b, DataType c, DataType d, DataType e) {

    auto b0 = 13.0 * (a - 2.0 * b + c) * (a - 2.0 * b + c) + 3.0 * (a - 4.0 * b + 3.0 * c) * (a - 4.0 * b + 3.0 * c);
    auto b1 = 13.0 * (b - 2.0 * c + d) * (b - 2.0 * c + d) + 3.0 * (b - d) * (b - d);
    auto b2 = 13.0 * (c - 2.0 * d + e) * (c - 2.0 * d + e) + 3.0 * (3.0 * c - 4.0 * d + e) * (3.0 * c - 4.0 * d + e);

    return weno_data{b0, b1, b2};
}

weno_data weno_solver::interpolation_p(DataType a, DataType b, DataType c, DataType d, DataType e) {

    auto f2 = (-a + 5.0 * b + 2.0 * c) / 6.0;
    auto f1 = (2.0 * b + 5.0 * c - d) / 6.0;
    auto f0 = (11.0 * c - 7.0 * d + 2.0 * e) / 6.0;
    return weno_data{f0, f1, f2};
}

weno_data weno_solver::interpolation_m(DataType a, DataType b, DataType c, DataType d, DataType e) {

    auto f2 = (-e + 5.0 * d + 2.0 * c) / 6.0;
    auto f1 = (2.0 * d + 5.0 * c - b) / 6.0;
    auto f0 = (11.0 * c - 7.0 * b + 2.0 * a) / 6.0;
    return weno_data{f0, f1, f2};
}

weno_data weno_solver::wenojs_ceofficients(DataType b0, DataType b1, DataType b2) {

// WENO-JS
//    auto a0 = 1.0 / (epsilon + b0) / (epsilon + b0);
//    auto a1 = 6.0 / (epsilon + b1) / (epsilon + b1);
//    auto a2 = 3.0 / (epsilon + b2) / (epsilon + b2);

// WENO-Z
    auto a0 = 1.0 * (1.0 + fabs(b0 - b2) / (epsilon + b0));
    auto a1 = 6.0 * (1.0 + fabs(b0 - b2) / (epsilon + b1));
    auto a2 = 3.0 * (1.0 + fabs(b0 - b2) / (epsilon + b2));

    auto w1 = a0 / (a0 + a1 + a2);
    auto w2 = a1 / (a0 + a1 + a2);
    auto w3 = a2 / (a0 + a1 + a2);
    return weno_data{w1, w2, w3};
}

void weno_solver::weno5_flux_x(scalar_data *f) {
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = -1; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);

                // left-biased stencil
                auto a = f->flux[index.i - 2][index.j][index.k];
                auto b = f->flux[index.i - 1][index.j][index.k];
                auto c = f->flux[index.i][index.j][index.k];
                auto d = f->flux[index.i + 1][index.j][index.k];
                auto e = f->flux[index.i + 2][index.j][index.k];

                auto fint = interpolation_m(a, b, c, d, e);
                auto indi_m = indicators_m(a, b, c, d, e);
                auto weights_m = wenojs_ceofficients(indi_m.x0, indi_m.x1, indi_m.x2);

                fm[index.i][index.j][index.k] =
                        weights_m.x0 * fint.x0 + weights_m.x1 * fint.x1 + weights_m.x2 * fint.x2;

                // right-biased stencil
                a = f->flux[index.i - 1][index.j][index.k];
                b = f->flux[index.i][index.j][index.k];
                c = f->flux[index.i + 1][index.j][index.k];
                d = f->flux[index.i + 2][index.j][index.k];
                e = f->flux[index.i + 3][index.j][index.k];

                fint = interpolation_p(a, b, c, d, e);
                auto indi_p = indicators_p(a, b, c, d, e);
                auto weights_p = wenojs_ceofficients(indi_p.x0, indi_p.x1, indi_p.x2);

                fp[index.i][index.j][index.k] =
                        weights_p.x0 * fint.x0 + weights_p.x1 * fint.x1 + weights_p.x2 * fint.x2;
            }
        }
    }
}

void weno_solver::weno5_flux_y(scalar_data *f) {
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = -1; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);

                // left-biased stencil
                auto a = f->flux[index.i][index.j - 2][index.k];
                auto b = f->flux[index.i][index.j - 1][index.k];
                auto c = f->flux[index.i][index.j][index.k];
                auto d = f->flux[index.i][index.j + 1][index.k];
                auto e = f->flux[index.i][index.j + 2][index.k];

                auto fint = interpolation_m(a, b, c, d, e);
                auto indi_m = indicators_m(a, b, c, d, e);
                auto weights_m = wenojs_ceofficients(indi_m.x0, indi_m.x1, indi_m.x2);

                fm[index.i][index.j][index.k] =
                        weights_m.x0 * fint.x0 + weights_m.x1 * fint.x1 + weights_m.x2 * fint.x2;

                // right-biased stencil
                a = f->flux[index.i][index.j - 1][index.k];
                b = f->flux[index.i][index.j][index.k];
                c = f->flux[index.i][index.j + 1][index.k];
                d = f->flux[index.i][index.j + 2][index.k];
                e = f->flux[index.i][index.j + 3][index.k];

                fint = interpolation_p(a, b, c, d, e);
                auto indi_p = indicators_p(a, b, c, d, e);
                auto weights_p = wenojs_ceofficients(indi_p.x0, indi_p.x1, indi_p.x2);

                fp[index.i][index.j][index.k] =
                        weights_p.x0 * fint.x0 + weights_p.x1 * fint.x1 + weights_p.x2 * fint.x2;
            }
        }
    }
}

void weno_solver::weno5_flux_z(scalar_data *f) {
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = -1; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);

                // left-biased stencil
                auto a = f->flux[index.i][index.j][index.k - 2];
                auto b = f->flux[index.i][index.j][index.k - 1];
                auto c = f->flux[index.i][index.j][index.k];
                auto d = f->flux[index.i][index.j][index.k + 1];
                auto e = f->flux[index.i][index.j][index.k + 2];

                auto fint = interpolation_m(a, b, c, d, e);
                auto indi_m = indicators_m(a, b, c, d, e);
                auto weights_m = wenojs_ceofficients(indi_m.x0, indi_m.x1, indi_m.x2);

                fm[index.i][index.j][index.k] =
                        weights_m.x0 * fint.x0 + weights_m.x1 * fint.x1 + weights_m.x2 * fint.x2;

                // right-biased stencil
                a = f->flux[index.i][index.j][index.k - 1];
                b = f->flux[index.i][index.j][index.k];
                c = f->flux[index.i][index.j][index.k + 1];
                d = f->flux[index.i][index.j][index.k + 2];
                e = f->flux[index.i][index.j][index.k + 3];

                fint = interpolation_p(a, b, c, d, e);
                auto indi_p = indicators_p(a, b, c, d, e);
                auto weights_p = wenojs_ceofficients(indi_p.x0, indi_p.x1, indi_p.x2);

                fp[index.i][index.j][index.k] =
                        weights_p.x0 * fint.x0 + weights_p.x1 * fint.x1 + weights_p.x2 * fint.x2;
            }
        }
    }
}

void weno_solver::weno5_find_fx(scalar_data *f, vector_data *vel) {
    weno5_flux_x(f);
#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = -1; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);
                if (vel->x.data[index.i][index.j][index.k] + vel->x.data[index.i + 1][index.j][index.k] > 0.0) {
                    fh[index.i][index.j][index.k] = fm[index.i][index.j][index.k];
                } else {
                    fh[index.i][index.j][index.k] = fp[index.i][index.j][index.k];
                }
            }
        }
    }
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);
                f->fx[index.i][index.j][index.k] =
                        (fh[index.i][index.j][index.k] - fh[index.i - 1][index.j][index.k]) / dx;
            }
        }
    }
}

void weno_solver::weno5_find_fy(scalar_data *f, vector_data *vel) {

    weno5_flux_y(f);
#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = -1; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);
                if (vel->y.data[index.i][index.j][index.k] + vel->y.data[index.i][index.j + 1][index.k] > 0.0) {
                    fh[index.i][index.j][index.k] = fm[index.i][index.j][index.k];
                } else {
                    fh[index.i][index.j][index.k] = fp[index.i][index.j][index.k];
                }
            }
        }
    }
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);
                f->fy[index.i][index.j][index.k] =
                        (fh[index.i][index.j][index.k] - fh[index.i][index.j - 1][index.k]) / dy;
            }
        }
    }

}

void weno_solver::weno5_find_fz(scalar_data *f, vector_data *vel) {

    weno5_flux_z(f);
#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = -1; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);
                if (vel->z.data[index.i][index.j][index.k] + vel->z.data[index.i][index.j][index.k + 1] > 0.0) {
                    fh[index.i][index.j][index.k] = fm[index.i][index.j][index.k];
                } else {
                    fh[index.i][index.j][index.k] = fp[index.i][index.j][index.k];
                }
            }
        }
    }
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i + 1, j + 1, k + 1);
                f->fz[index.i][index.j][index.k] =
                        (fh[index.i][index.j][index.k] - fh[index.i][index.j][index.k - 1]) / dz;
            }
        }
    }

}

void weno_solver::weno5_find_derivatives(scalar_data *f, vector_data *vel) {
    weno5_find_fx(f, vel);
    if (f->ndim > 1) {
        weno5_find_fy(f, vel);
    }
    if (f->ndim > 2) {
        weno5_find_fz(f, vel);
    }
}
