//
// Created by leo on 10/4/22.
//

#include "wrapper_func.h"

void find_heavyside(wrapper *lsf) {
#pragma omp parallel for default(none) shared(lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx; ++i) {
        for (int j = 0; j < lsf->scalar->Ny; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {
                lsf->dummy->heaviside[i][j][k] = Heaviside(lsf->scalar->data[i][j][k], lsf->params->ls_width);
            }
        }
    }
}

void find_delta(wrapper *lsf) {
#pragma omp parallel for default(none) shared(lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx; ++i) {
        for (int j = 0; j < lsf->scalar->Ny; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {
                lsf->dummy->delta[i][j][k] = Delta(lsf->scalar->data[i][j][k], lsf->params->ls_width);
            }
        }
    }
}

void find_sign(wrapper *lsf) {
#pragma omp parallel for default(none) shared(lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx; ++i) {
        for (int j = 0; j < lsf->scalar->Ny; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {
                lsf->dummy->sign[i][j][k] = Sign(lsf->scalar->data[i][j][k], lsf->params->ls_width);
            }
        }
    }
}

void find_gradient(wrapper *lsf) {
#pragma omp parallel for default(none) shared(lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx; ++i) {
        for (int j = 0; j < lsf->scalar->Ny; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {
                auto gradient = pow(lsf->scalar->fx[i][j][k], 2);
                if (lsf->scalar->ndim > 1) {
                    gradient += pow(lsf->scalar->fy[i][j][k], 2);
                }
                if (lsf->scalar->ndim > 2) {
                    gradient += pow(lsf->scalar->fz[i][j][k], 2);
                }
                lsf->dummy->grad[i][j][k] = sqrt(gradient);
            }
        }
    }
}

void store(wrapper *f) {
    if (f->is_scalar) {
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    f->dummy->tmp[i][j][k] = f->scalar->data[i][j][k];
                }
            }
        }
    } else {
        for (int i = 0; i < f->vector->x.Nx; ++i) {
            for (int j = 0; j < f->vector->x.Ny; ++j) {
                for (int k = 0; k < f->vector->x.Nz; ++k) {
                    f->dummy->u_tmp[i][j][k] = f->vector->x.data[i][j][k];
                    f->dummy->v_tmp[i][j][k] = f->vector->y.data[i][j][k];
                    f->dummy->w_tmp[i][j][k] = f->vector->z.data[i][j][k];
                }
            }
        }
    }
}

DataType l2norm(wrapper *f) {

    DataType error = 0.0;

    if (f->is_scalar) {
#pragma omp parallel for default(none) shared(f) reduction(+:error) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    error += pow(f->dummy->tmp[i][j][k] - f->scalar->data[i][j][k], 2);
                }
            }
        }
        error = sqrt(error / (f->scalar->Nx * f->scalar->Ny * f->scalar->Nz));
    } else {
#pragma omp parallel for default(none) shared(f) reduction(+:error) collapse(3)
        for (int i = 0; i < f->vector->x.Nx; ++i) {
            for (int j = 0; j < f->vector->x.Ny; ++j) {
                for (int k = 0; k < f->vector->x.Nz; ++k) {
                    error += pow(f->dummy->u_tmp[i][j][k] - f->vector->x.data[i][j][k], 2);
                    error += pow(f->dummy->v_tmp[i][j][k] - f->vector->y.data[i][j][k], 2);
                    error += pow(f->dummy->w_tmp[i][j][k] - f->vector->z.data[i][j][k], 2);
                }
            }
        }
        error = sqrt(error / (f->vector->x.Nx * f->vector->x.Ny * f->vector->x.Nz));
    }

    return error;
}

void integrate_a(wrapper *f) {

    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->dummy->a_int[i][j][k] = 0.0;
            }
        }
    }

    if (f->scalar->ndim == 1) {
#pragma omp parallel for default(none) shared(f) collapse(3)
        for (int i = 0; i < f->scalar->nx; ++i) {
            for (int j = 0; j < f->scalar->ny; ++j) {
                for (int k = 0; k < f->scalar->nz; ++k) {
                    auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                    // Simplson's rule
                    auto f0 =
                            0.5 * (f->dummy->a[index.i][index.j][index.k] + f->dummy->a[index.i - 1][index.j][index.k]);
                    auto f1 = f->dummy->a[index.i][index.j][index.k];
                    auto f2 =
                            0.5 * (f->dummy->a[index.i][index.j][index.k] + f->dummy->a[index.i + 1][index.j][index.k]);
                    f->dummy->a_int[index.i][index.j][index.k] = (f0 + 4.0 * f1 + f2) / 6.0;
                }
            }
        }
    } else if (f->scalar->ndim == 2) {
#pragma omp parallel for default(none) shared(f) collapse(3)
        for (int i = 0; i < f->scalar->nx; ++i) {
            for (int j = 0; j < f->scalar->ny; ++j) {
                for (int k = 0; k < f->scalar->nz; ++k) {
                    auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                    // Simplson's rule
                    f->dummy->a_int[index.i][index.j][index.k] = 16.0 * f->dummy->a[index.i][index.j][index.k];
                    for (int ii = 0; ii < 2; ++ii) {
                        for (int jj = 0; jj < 2; ++jj) {
                            f->dummy->a_int[index.i][index.j][index.k] +=
                                    f->dummy->a[index.i - 1 + 2 * ii][index.j - 1 + 2 * jj][index.k];
                        }
                    }
                    f->dummy->a_int[index.i][index.j][index.k] /= 24.0;
                }
            }
        }
    } else {
        for (int i = 0; i < f->scalar->nx; ++i) {
            for (int j = 0; j < f->scalar->ny; ++j) {
                for (int k = 0; k < f->scalar->nz; ++k) {
                    auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                    // Simplson's rule
                    f->dummy->a_int[index.i][index.j][index.k] = 51.0 * f->dummy->a[index.i][index.j][index.k];
                    for (int ii = 0; ii < 2; ++ii) {
                        for (int jj = 0; jj < 2; ++jj) {
                            for (int kk = 0; kk < 2; ++kk) {
                                f->dummy->a_int[index.i][index.j][index.k] +=
                                        f->dummy->a[index.i - 1 + 2 * ii][index.j - 1 + 2 * jj][index.k - 1 + 2 * kk];
                            }
                        }
                    }
                    f->dummy->a_int[index.i][index.j][index.k] /= 77.0;
                }
            }
        }
    }
}

void integrate_b(wrapper *f) {

    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                f->dummy->b_int[i][j][k] = 0.0;
            }
        }
    }

    if (f->scalar->ndim == 1) {
#pragma omp parallel for default(none) shared(f) collapse(3)
        for (int i = 0; i < f->scalar->nx; ++i) {
            for (int j = 0; j < f->scalar->ny; ++j) {
                for (int k = 0; k < f->scalar->nz; ++k) {
                    auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                    // Simplson's rule
                    auto f0 =
                            0.5 * (f->dummy->b[index.i][index.j][index.k] + f->dummy->b[index.i - 1][index.j][index.k]);
                    auto f1 = f->dummy->b[index.i][index.j][index.k];
                    auto f2 =
                            0.5 * (f->dummy->b[index.i][index.j][index.k] + f->dummy->b[index.i + 1][index.j][index.k]);
                    f->dummy->b_int[index.i][index.j][index.k] = (f0 + 4.0 * f1 + f2) / 6.0;
                }
            }
        }
    } else if (f->scalar->ndim == 2) {
#pragma omp parallel for default(none) shared(f) collapse(3)
        for (int i = 0; i < f->scalar->nx; ++i) {
            for (int j = 0; j < f->scalar->ny; ++j) {
                for (int k = 0; k < f->scalar->nz; ++k) {
                    auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                    // Simplson's rule
                    f->dummy->b_int[index.i][index.j][index.k] = 16.0 * f->dummy->b[index.i][index.j][index.k];
                    for (int ii = 0; ii < 2; ++ii) {
                        for (int jj = 0; jj < 2; ++jj) {
                            f->dummy->b_int[index.i][index.j][index.k] +=
                                    f->dummy->b[index.i - 1 + 2 * ii][index.j - 1 + 2 * jj][index.k];
                        }
                    }
                    f->dummy->b_int[index.i][index.j][index.k] /= 24.0;
                }
            }
        }
    } else {
        for (int i = 0; i < f->scalar->nx; ++i) {
            for (int j = 0; j < f->scalar->ny; ++j) {
                for (int k = 0; k < f->scalar->nz; ++k) {
                    auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                    // Simplson's rule
                    f->dummy->b_int[index.i][index.j][index.k] = 51.0 * f->dummy->b[index.i][index.j][index.k];
                    for (int ii = 0; ii < 2; ++ii) {
                        for (int jj = 0; jj < 2; ++jj) {
                            for (int kk = 0; kk < 2; ++kk) {
                                f->dummy->b_int[index.i][index.j][index.k] +=
                                        f->dummy->b[index.i - 1 + 2 * ii][index.j - 1 + 2 * jj][index.k - 1 + 2 * kk];
                            }
                        }
                    }
                }
            }
        }
    }
}
