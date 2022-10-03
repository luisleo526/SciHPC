//
// Created by leo on 10/4/22.
//

#include "godunov_gradient.h"

DataType godunov_limiter_p(DataType fp, DataType fm) {
    return std::pow(std::max(-std::min(fp, 0.0), std::max(fm, 0.0)), 2);
}

DataType godunov_limiter_m(DataType fp, DataType fm) {
    return std::pow(std::max(std::max(fp, 0.0), -std::min(fm, 0.0)), 2);
}

void godunov_gradient(wrapper *f, structured_grid *geo) {

    // Initialized for gradient
#pragma omp parallel for default(none) shared(f, geo) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->dummy->grad[i][j][k] = 0.0;
            }
        }
    }

    // x direction
    // prepare first order derivative
#pragma omp parallel for default(none) shared(f, geo) collapse(3)
    for (int i = 1; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->scalar->flux[i][j][k] = (f->scalar->data[i][j][k] - f->scalar->data[i - 1][j][k]) / geo->dx;
            }
        }
    }
    // boundary condition for first order derivative
#pragma omp parallel for default(none) shared(f, geo) collapse(2)
    for (int j = 0; j < f->scalar->Ny; ++j) {
        for (int k = 0; k < f->scalar->Nz; ++k) {
            auto index_l = f->scalar->index_mapping(1, 1, 1);
            auto index_r = f->scalar->index_mapping(f->scalar->nx, 1, 1);
            for (int i = 1; i <= f->scalar->ghc; ++i) {
                f->scalar->flux[index_l.i - i][j][k] = f->scalar->flux[index_l.i][j][k];
                f->scalar->flux[index_r.i + i][j][k] = f->scalar->flux[index_r.i][j][k];
            }
        }
    }
    // find fluxes in x direction
    f->solvers->weno->wenojs_flux_x(f->scalar);
#pragma omp parallel for default(none) shared(f, geo) collapse(3)
    for (int i = 0; i < f->scalar->nx; ++i) {
        for (int j = 0; j < f->scalar->ny; ++j) {
            for (int k = 0; k < f->scalar->nz; ++k) {
                auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                if (f->dummy->sign[index.i][index.j][index.k] > 0.0) {
                    f->dummy->grad[index.i][index.j][index.k] = godunov_limiter_p(
                            f->solvers->weno->fp[index.i][index.j][index.k],
                            f->solvers->weno->fm[index.i][index.j][index.k]);
                } else {
                    f->dummy->grad[index.i][index.j][index.k] = godunov_limiter_m(
                            f->solvers->weno->fp[index.i][index.j][index.k],
                            f->solvers->weno->fm[index.i][index.j][index.k]);
                }
            }
        }
    }

    if (f->scalar->ndim > 1) {
        // y direction
        // prepare first order derivative
#pragma omp parallel for default(none) shared(f, geo) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 1; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    f->scalar->flux[i][j][k] = (f->scalar->data[i][j][k] - f->scalar->data[i][j - 1][k]) / geo->dy;
                }
            }
        }

        // boundary condition for first order derivative
#pragma omp parallel for default(none) shared(f, geo) collapse(2)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                auto index_l = f->scalar->index_mapping(1, 1, 1);
                auto index_r = f->scalar->index_mapping(1, f->scalar->ny, 1);
                for (int j = 1; j <= f->scalar->ghc; ++j) {
                    f->scalar->flux[i][index_l.j - j][k] = f->scalar->flux[i][index_l.j][k];
                    f->scalar->flux[i][index_r.j + j][k] = f->scalar->flux[i][index_r.j][k];
                }
            }
        }

        // find fluxes in y direction
        f->solvers->weno->wenojs_flux_y(f->scalar);
#pragma omp parallel for default(none) shared(f, geo) collapse(3)
        for (int i = 0; i < f->scalar->nx; ++i) {
            for (int j = 0; j < f->scalar->ny; ++j) {
                for (int k = 0; k < f->scalar->nz; ++k) {
                    auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                    if (f->dummy->sign[index.i][index.j][index.k] > 0.0) {
                        f->dummy->grad[index.i][index.j][index.k] += godunov_limiter_p(
                                f->solvers->weno->fp[index.i][index.j][index.k],
                                f->solvers->weno->fm[index.i][index.j][index.k]);
                    } else {
                        f->dummy->grad[index.i][index.j][index.k] += godunov_limiter_m(
                                f->solvers->weno->fp[index.i][index.j][index.k],
                                f->solvers->weno->fm[index.i][index.j][index.k]);
                    }
                }
            }
        }
    }

    if (f->scalar->ndim > 2) {
        // z direction
        // prepare first order derivative
#pragma omp parallel for default(none) shared(f, geo) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 1; k < f->scalar->Nz; ++k) {
                    f->scalar->flux[i][j][k] = (f->scalar->data[i][j][k] - f->scalar->data[i][j][k - 1]) / geo->dz;
                }
            }
        }

        // boundary condition for first order derivative
#pragma omp parallel for default(none) shared(f, geo) collapse(2)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                auto index_l = f->scalar->index_mapping(1, 1, 1);
                auto index_r = f->scalar->index_mapping(1, 1, f->scalar->nz);
                for (int k = 1; k <= f->scalar->ghc; ++k) {
                    f->scalar->flux[i][j][index_l.k - k] = f->scalar->flux[i][j][index_l.k];
                    f->scalar->flux[i][j][index_r.k + k] = f->scalar->flux[i][j][index_r.k];
                }
            }
        }

        // find fluxes in z direction
        f->solvers->weno->wenojs_flux_z(f->scalar);
#pragma omp parallel for default(none) shared(f, geo) collapse(3)
        for (int i = 0; i < f->scalar->nx; ++i) {
            for (int j = 0; j < f->scalar->ny; ++j) {
                for (int k = 0; k < f->scalar->nz; ++k) {
                    auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                    if (f->dummy->sign[index.i][index.j][index.k] > 0.0) {
                        f->dummy->grad[index.i][index.j][index.k] += godunov_limiter_p(
                                f->solvers->weno->fp[index.i][index.j][index.k],
                                f->solvers->weno->fm[index.i][index.j][index.k]);
                    } else {
                        f->dummy->grad[index.i][index.j][index.k] += godunov_limiter_m(
                                f->solvers->weno->fp[index.i][index.j][index.k],
                                f->solvers->weno->fm[index.i][index.j][index.k]);
                    }
                }
            }
        }
    }

    // finalize for gradient
#pragma omp parallel for default(none) shared(f, geo) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->dummy->grad[i][j][k] = sqrt(f->dummy->grad[i][j][k]);
            }
        }
    }
}

void stabilized_upon_gradient(wrapper *f, structured_grid *geo) {
    godunov_gradient(f, geo);
    DataType max_grad = 0.0;

#pragma omp parallel for reduction(max:max_grad) default(none) shared(f) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                max_grad = std::max(max_grad, f->dummy->grad[i][j][k]);
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



