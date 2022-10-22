//
// Created by leo on 10/4/22.
//

#include "godunov_gradient.h"

void godunov_gradient(wrapper *f) {

    // Initialized for gradient
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->dummy->grad[i][j][k] = 0.0;
            }
        }
    }

    // x direction

    auto phi = f->scalar->data;

#pragma omp parallel for default(none) shared(f, phi) collapse(3)
    for (int I = 0; I < f->scalar->nx; ++I) {
        for (int J = 0; J < f->scalar->ny; ++J) {
            for (int K = 0; K < f->scalar->nz; ++K) {
                auto index = f->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                auto v = 1.0 / (12.0 * f->geo->dx) * (-(phi[i - 1][j][k] - phi[i - 2][j][k])
                                                      + 7.0 * (phi[i][j][k] - phi[i - 1][j][k])
                                                      + 7.0 * (phi[i + 1][j][k] - phi[i][j][k])
                                                      - (phi[i + 2][j][k] - phi[i + 1][j][k]));

                auto up = v + 1.0 / f->geo->dx *
                              weno5_for_godunov(phi[i + 3][j][k] - 2.0 * phi[i + 2][j][k] + phi[i + 1][j][k],
                                                phi[i + 2][j][k] - 2.0 * phi[i + 1][j][k] + phi[i][j][k],
                                                phi[i + 1][j][k] - 2.0 * phi[i][j][k] + phi[i - 1][j][k],
                                                phi[i][j][k] - 2.0 * phi[i - 1][j][k] + phi[i - 2][j][k]);

                auto um = v - 1.0 / f->geo->dx *
                              weno5_for_godunov(phi[i - 3][j][k] - 2.0 * phi[i - 2][j][k] + phi[i - 1][j][k],
                                                phi[i - 2][j][k] - 2.0 * phi[i - 1][j][k] + phi[i][j][k],
                                                phi[i - 1][j][k] - 2.0 * phi[i][j][k] + phi[i + 1][j][k],
                                                phi[i][j][k] - 2.0 * phi[i + 1][j][k] + phi[i + 2][j][k]);

                if (f->dummy->sign[i][j][k] > 0.0) {
                    auto upm = -fmin(up, 0.0);
                    auto ump = fmax(um, 0.0);
                    f->dummy->grad[i][j][k] = pow(fmax(upm, ump), 2);
                } else {
                    auto upp = fmax(up, 0.0);
                    auto upm = -fmin(um, 0.0);
                    f->dummy->grad[i][j][k] = pow(fmax(upp, upm), 2);
                }

            }
        }
    }

    if (f->scalar->ndim > 1) {
        // y direction
#pragma omp parallel for default(none) shared(f, phi) collapse(3)
        for (int I = 0; I < f->scalar->nx; ++I) {
            for (int J = 0; J < f->scalar->ny; ++J) {
                for (int K = 0; K < f->scalar->nz; ++K) {
                    auto index = f->scalar->index_mapping(I + 1, J + 1, K + 1);
                    auto i = index.i;
                    auto j = index.j;
                    auto k = index.k;

                    auto v = 1.0 / (12.0 * f->geo->dy) * (-(phi[i][j - 1][k] - phi[i][j - 2][k])
                                                          + 7.0 * (phi[i][j][k] - phi[i][j - 1][k])
                                                          + 7.0 * (phi[i][j + 1][k] - phi[i][j][k])
                                                          - (phi[i][j + 2][k] - phi[i][j + 1][k]));

                    auto vp = v + 1.0 / f->geo->dy *
                                  weno5_for_godunov(phi[i][j + 3][k] - 2.0 * phi[i][j + 2][k] + phi[i][j + 1][k],
                                                    phi[i][j + 2][k] - 2.0 * phi[i][j + 1][k] + phi[i][j][k],
                                                    phi[i][j + 1][k] - 2.0 * phi[i][j][k] + phi[i][j - 1][k],
                                                    phi[i][j][k] - 2.0 * phi[i][j - 1][k] + phi[i][j - 2][k]);

                    auto vm = v - 1.0 / f->geo->dy *
                                  weno5_for_godunov(phi[i][j - 3][k] - 2.0 * phi[i][j - 2][k] + phi[i][j - 1][k],
                                                    phi[i][j - 2][k] - 2.0 * phi[i][j - 1][k] + phi[i][j][k],
                                                    phi[i][j - 1][k] - 2.0 * phi[i][j][k] + phi[i][j + 1][k],
                                                    phi[i][j][k] - 2.0 * phi[i][j + 1][k] + phi[i][j + 2][k]);

                    if (f->dummy->sign[i][j][k] > 0.0) {
                        auto vpm = -fmin(vp, 0.0);
                        auto vmp = fmax(vm, 0.0);
                        f->dummy->grad[i][j][k] += pow(fmax(vpm, vmp), 2);
                    } else {
                        auto vpp = fmax(vp, 0.0);
                        auto vpm = -fmin(vm, 0.0);
                        f->dummy->grad[i][j][k] += pow(fmax(vpp, vpm), 2);
                    }
                }
            }
        }
    }

    if (f->scalar->ndim > 2) {
        // z direction
#pragma omp parallel for default(none) shared(f, phi) collapse(3)
        for (int I = 0; I < f->scalar->nx; ++I) {
            for (int J = 0; J < f->scalar->ny; ++J) {
                for (int K = 0; K < f->scalar->nz; ++K) {
                    auto index = f->scalar->index_mapping(I + 1, J + 1, K + 1);
                    auto i = index.i;
                    auto j = index.j;
                    auto k = index.k;

                    auto v = 1.0 / (12.0 * f->geo->dz) * (-(phi[i][j][k - 1] - phi[i][j][k - 2])
                                                          + 7.0 * (phi[i][j][k] - phi[i][j][k - 1])
                                                          + 7.0 * (phi[i][j][k + 1] - phi[i][j][k])
                                                          - (phi[i][j][k + 2] - phi[i][j][k + 1]));

                    auto wp = v + 1.0 / f->geo->dz *
                                  weno5_for_godunov(phi[i][j][k + 3] - 2.0 * phi[i][j][k + 2] + phi[i][j][k + 1],
                                                    phi[i][j][k + 2] - 2.0 * phi[i][j][k + 1] + phi[i][j][k],
                                                    phi[i][j][k + 1] - 2.0 * phi[i][j][k] + phi[i][j][k - 1],
                                                    phi[i][j][k] - 2.0 * phi[i][j][k - 1] + phi[i][j][k - 2]);

                    auto wm = v - 1.0 / f->geo->dz *
                                  weno5_for_godunov(phi[i][j][k - 3] - 2.0 * phi[i][j][k - 2] + phi[i][j][k - 1],
                                                    phi[i][j][k - 2] - 2.0 * phi[i][j][k - 1] + phi[i][j][k],
                                                    phi[i][j][k - 1] - 2.0 * phi[i][j][k] + phi[i][j][k + 1],
                                                    phi[i][j][k] - 2.0 * phi[i][j][k + 1] + phi[i][j][k + 2]);

                    if (f->dummy->sign[i][j][k] > 0.0) {
                        auto wpm = -fmin(wp, 0.0);
                        auto wmp = fmax(wm, 0.0);
                        f->dummy->grad[i][j][k] += pow(fmax(wpm, wmp), 2);
                    } else {
                        auto wpp = fmax(wp, 0.0);
                        auto wpm = -fmin(wm, 0.0);
                        f->dummy->grad[i][j][k] += pow(fmax(wpp, wpm), 2);
                    }

                }
            }
        }
    }

    // finalize for gradient
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->dummy->grad[i][j][k] = sqrt(f->dummy->grad[i][j][k]);
            }
        }
    }
}

void stabilized_upon_gradient(wrapper *f) {
    godunov_gradient(f);
    DataType max_grad = 0.0;

#pragma omp parallel for reduction(max:max_grad) default(none) shared(f) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                max_grad = fmax(max_grad, f->dummy->grad[i][j][k]);
            }
        }
    }

#pragma omp parallel for default(none) shared(f, max_grad) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->scalar->data[i][j][k] /= max_grad;
            }
        }
    }

}

DataType weno5_for_godunov(DataType a, DataType b, DataType c, DataType d) {
    auto eps = epsilon;
    auto is0 = 13.0 * pow(a - b, 2.0) + 3.0 * pow(a - 3.0 * b, 2.0);
    auto is1 = 13.0 * pow(b - c, 2.0) + 3.0 * pow(b + c, 2.0);
    auto is2 = 13.0 * pow(c - d, 2.0) + 3.0 * pow(3.0 * c - d, 2.0);
    auto alp0 = 1.0 / pow(eps + is0, 2.0);
    auto alp1 = 6.0 / pow(eps + is1, 2.0);
    auto alp2 = 3.0 / pow(eps + is2, 2.0);
//    auto alp0 = 1.0 * (1.0 + fabs(is0 - is2) / (epsilon + is0));
//    auto alp1 = 6.0 * (1.0 + fabs(is0 - is2) / (epsilon + is1));
//    auto alp2 = 3.0 * (1.0 + fabs(is0 - is2) / (epsilon + is2));

    auto w0 = alp0 / (alp0 + alp1 + alp2);
    auto w2 = alp2 / (alp0 + alp1 + alp2);
    return w0 / 3.0 * (a - 2.0 * b + c) + (w2 - 0.5) / 6.0 * (b - 2.0 * c + d);
}



