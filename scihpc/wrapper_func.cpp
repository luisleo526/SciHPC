//
// Created by leo on 10/4/22.
//

#include "wrapper_func.h"
#include <iostream>

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

void store_tmp(wrapper *f) {
    if (f->is_scalar) {
#pragma omp parallel for default(none) shared(f) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    f->dummy->tmp[i][j][k] = f->scalar->data[i][j][k];
                }
            }
        }
    } else {
#pragma omp parallel for default(none) shared(f) collapse(3)
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

void find_curvature(wrapper *lsf) {

    if (lsf->scalar->ndim == 2) {
#pragma omp parallel for default(none) shared(lsf) collapse(3)
        for (int i = 0; i < lsf->scalar->Nx; ++i) {
            for (int j = 0; j < lsf->scalar->Ny; ++j) {
                for (int k = 0; k < lsf->scalar->Nz; ++k) {
                    auto fx = lsf->scalar->fx[i][j][k];
                    auto fy = lsf->scalar->fy[i][j][k];
                    auto fxx = lsf->scalar->fxx[i][j][k];
                    auto fyy = lsf->scalar->fyy[i][j][k];
                    auto fxy = lsf->scalar->fxy[i][j][k];
                    lsf->dummy->curvature[i][j][k] = -(fxx * fy * fy + fyy * fx * fx - 2.0 * fxy * fx * fy) /
                                                     pow(fx * fx + fy * fy, 1.5);

                    //Correction
                    lsf->dummy->curvature[i][j][k] =
                            1.0 / (1.0 / lsf->dummy->curvature[i][j][k] + lsf->scalar->data[i][j][k]);
                }
            }
        }
    } else if (lsf->scalar->ndim == 3) {
#pragma omp parallel for default(none) shared(lsf) collapse(3)
        for (int i = 0; i < lsf->scalar->Nx; ++i) {
            for (int j = 0; j < lsf->scalar->Ny; ++j) {
                for (int k = 0; k < lsf->scalar->Nz; ++k) {
                    auto fx = lsf->scalar->fx[i][j][k];
                    auto fy = lsf->scalar->fy[i][j][k];
                    auto fz = lsf->scalar->fz[i][j][k];
                    auto fxx = lsf->scalar->fxx[i][j][k];
                    auto fyy = lsf->scalar->fyy[i][j][k];
                    auto fzz = lsf->scalar->fzz[i][j][k];
                    auto fxy = lsf->scalar->fxy[i][j][k];
                    auto fxz = lsf->scalar->fzx[i][j][k];
                    auto fyz = lsf->scalar->fyz[i][j][k];
                    lsf->dummy->curvature[i][j][k] =
                            -(fx * fx * (fyy + fzz) + fy * fy * (fxx + fzz) + fz * fz * (fxx + fyy)
                              - 2.0 * (fx * fy * fxy + fx * fz * fxz + fy * fz * fyz)) /
                            pow(fx * fx + fy * fy + fz * fz, 1.5);

                    //Correction
                    lsf->dummy->curvature[i][j][k] =
                            1.0 / (2.0 / lsf->dummy->curvature[i][j][k] + lsf->scalar->data[i][j][k]);
                }
            }
        }
    }
}

void all_to_face_x(wrapper *ref, wrapper *tgt) {
    if (ref->vector->x.ndim == 1) {
#pragma omp parallel for default(none) shared(ref, tgt) collapse(3)
        for (int i = 1; i < ref->vector->x.Nx - 1; ++i) {
            for (int j = 0; j < ref->vector->x.Ny; ++j) {
                for (int k = 0; k < ref->vector->x.Nz; ++k) {
                    tgt->vector->x.data[i][j][k] = ref->vector->x.data[i][j][k];
                }
            }
        }
    } else if (ref->vector->x.ndim == 2) {
#pragma omp parallel for default(none) shared(ref, tgt) collapse(3)
        for (int i = 1; i < ref->vector->x.Nx - 1; ++i) {
            for (int j = 1; j < ref->vector->x.Ny - 1; ++j) {
                for (int k = 0; k < ref->vector->x.Nz; ++k) {
                    tgt->vector->x.data[i][j][k] = ref->vector->x.data[i][j][k];
                    tgt->vector->y.data[i][j][k] =
                            0.25 * (ref->vector->y.data[i][j][k] + ref->vector->y.data[i + 1][j][k] +
                                    ref->vector->y.data[i][j - 1][k] + ref->vector->y.data[i + 1][j - 1][k]);
                }
            }
        }
    } else if (ref->vector->x.ndim == 3) {
#pragma omp parallel for default(none) shared(ref, tgt) collapse(3)
        for (int i = 1; i < ref->vector->x.Nx - 1; ++i) {
            for (int j = 1; j < ref->vector->x.Ny - 1; ++j) {
                for (int k = 1; k < ref->vector->x.Nz - 1; ++k) {
                    tgt->vector->x.data[i][j][k] = ref->vector->x.data[i][j][k];
                    tgt->vector->y.data[i][j][k] =
                            0.25 * (ref->vector->y.data[i][j][k] + ref->vector->y.data[i + 1][j][k] +
                                    ref->vector->y.data[i][j - 1][k] + ref->vector->y.data[i + 1][j - 1][k]);
                    tgt->vector->z.data[i][j][k] =
                            0.25 * (ref->vector->z.data[i][j][k] + ref->vector->z.data[i + 1][j][k] +
                                    ref->vector->z.data[i][j][k - 1] + ref->vector->z.data[i + 1][j][k - 1]);
                }
            }
        }
    }
}

void all_to_face_y(wrapper *ref, wrapper *tgt) {
    if (ref->vector->y.ndim == 2) {
#pragma omp parallel for default(none) shared(ref, tgt) collapse(3)
        for (int i = 1; i < ref->vector->y.Nx - 1; ++i) {
            for (int j = 1; j < ref->vector->y.Ny - 1; ++j) {
                for (int k = 0; k < ref->vector->y.Nz; ++k) {
                    tgt->vector->y.data[i][j][k] = ref->vector->y.data[i][j][k];
                    tgt->vector->x.data[i][j][k] =
                            0.25 * (ref->vector->x.data[i][j][k] + ref->vector->x.data[i][j + 1][k] +
                                    ref->vector->x.data[i - 1][j][k] + ref->vector->x.data[i - 1][j + 1][k]);
                }
            }
        }
    } else if (ref->vector->y.ndim == 3) {
#pragma omp parallel for default(none) shared(ref, tgt) collapse(3)
        for (int i = 1; i < ref->vector->y.Nx - 1; ++i) {
            for (int j = 1; j < ref->vector->y.Ny - 1; ++j) {
                for (int k = 1; k < ref->vector->y.Nz - 1; ++k) {
                    tgt->vector->y.data[i][j][k] = ref->vector->y.data[i][j][k];
                    tgt->vector->x.data[i][j][k] =
                            0.25 * (ref->vector->x.data[i][j][k] + ref->vector->x.data[i][j + 1][k] +
                                    ref->vector->x.data[i - 1][j][k] + ref->vector->x.data[i - 1][j + 1][k]);
                    tgt->vector->z.data[i][j][k] =
                            0.25 * (ref->vector->z.data[i][j][k] + ref->vector->z.data[i][j + 1][k] +
                                    ref->vector->z.data[i][j][k - 1] + ref->vector->z.data[i][j + 1][k - 1]);
                }
            }
        }
    }
}

void all_to_face_z(wrapper *ref, wrapper *tgt) {
    if (ref->vector->z.ndim == 3) {
#pragma omp parallel for default(none) shared(ref, tgt) collapse(3)
        for (int i = 1; i < ref->vector->z.Nx - 1; ++i) {
            for (int j = 1; j < ref->vector->z.Ny - 1; ++j) {
                for (int k = 1; k < ref->vector->z.Nz - 1; ++k) {
                    tgt->vector->z.data[i][j][k] = ref->vector->z.data[i][j][k];
                    tgt->vector->x.data[i][j][k] =
                            0.25 * (ref->vector->x.data[i][j][k] + ref->vector->x.data[i][j][k + 1] +
                                    ref->vector->x.data[i - 1][j][k] + ref->vector->x.data[i - 1][j][k + 1]);
                    tgt->vector->y.data[i][j][k] =
                            0.25 * (ref->vector->y.data[i][j][k] + ref->vector->y.data[i][j][k + 1] +
                                    ref->vector->y.data[i][j - 1][k] + ref->vector->y.data[i][j - 1][k + 1]);
                }
            }
        }
    }
}

void node_from_face(wrapper *ref, wrapper *tgt) {
    if (ref->vector->x.ndim == 2) {
        for (int i = 1; i < ref->vector->x.Nx - 1; ++i) {
            for (int j = 1; j < ref->vector->x.Ny - 1; ++j) {
                for (int k = 0; k < ref->vector->x.Nz; ++k) {
                    tgt->vector->x.data[i][j][k] =
                            0.5 * (ref->vector->x.data[i][j][k] + ref->vector->x.data[i - 1][j][k]);
                    tgt->vector->y.data[i][j][k] =
                            0.5 * (ref->vector->y.data[i][j][k] + ref->vector->y.data[i][j - 1][k]);
                }
            }
        }
    } else {
        for (int i = 1; i < ref->vector->x.Nx - 1; ++i) {
            for (int j = 1; j < ref->vector->x.Ny - 1; ++j) {
                for (int k = 1; k < ref->vector->x.Nz; ++k) {
                    tgt->vector->x.data[i][j][k] =
                            0.5 * (ref->vector->x.data[i][j][k] + ref->vector->x.data[i - 1][j][k]);
                    tgt->vector->y.data[i][j][k] =
                            0.5 * (ref->vector->y.data[i][j][k] + ref->vector->y.data[i][j - 1][k]);
                    tgt->vector->z.data[i][j][k] =
                            0.5 * (ref->vector->z.data[i][j][k] + ref->vector->z.data[i][j][k - 1]);
                }
            }
        }
    }
}

DataType divergence(wrapper *vel) {
    DataType max_div;

    if (vel->vector->x.ndim == 2) {
        max_div = 0.0;
#pragma omp parallel for default(none) shared(vel) reduction(max:max_div) collapse(3)
        for (int i = 0; i < vel->vector->x.nx; ++i) {
            for (int j = 0; j < vel->vector->x.ny; ++j) {
                for (int k = 0; k < vel->vector->x.Nz; ++k) {
                    auto index = vel->vector->x.index_mapping(i + 1, j + 1, k + 1);
                    auto I = index.i;
                    auto J = index.j;
                    auto K = index.k;
                    auto div = (vel->vector->x.data[I][J][K] - vel->vector->x.data[I - 1][J][K]) / vel->geo->dx +
                               (vel->vector->y.data[I][J][K] - vel->vector->y.data[I][J - 1][K]) / vel->geo->dy;
                    max_div = fmax(max_div, fabs(div));
                }
            }
        }
    } else {
        max_div = 0.0;
#pragma omp parallel for default(none) shared(vel) reduction(max:max_div) collapse(3)
        for (int i = 0; i < vel->vector->x.nx; ++i) {
            for (int j = 0; j < vel->vector->x.ny; ++j) {
                for (int k = 0; k < vel->vector->x.nz; ++k) {
                    auto index = vel->vector->x.index_mapping(i + 1, j + 1, k + 1);
                    auto I = index.i;
                    auto J = index.j;
                    auto K = index.k;
                    auto div = (vel->vector->x.data[I][J][K] - vel->vector->x.data[I - 1][J][K]) / vel->geo->dx +
                               (vel->vector->y.data[I][J][K] - vel->vector->y.data[I][J - 1][K]) / vel->geo->dy +
                               (vel->vector->z.data[I][J][K] - vel->vector->z.data[I][J][K - 1]) / vel->geo->dz;
                    max_div = fmax(max_div, fabs(div));
                }
            }
        }
    }
    return max_div;
}

void find_dt(wrapper *vel) {

    DataType max_u, max_v, max_w, max_curvature;
    max_u = max_v = max_w = 0.0;
#pragma omp parallel for default(none) shared(vel) reduction(max:max_u, max_v, max_w, max_curvature) collapse(3)
    for (int I = 0; I < vel->vector->x.nx; ++I) {
        for (int J = 0; J < vel->vector->x.ny; ++J) {
            for (int K = 0; K < vel->vector->x.nz; ++K) {
                auto index = vel->vector->x.index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;
                max_u = fmax(max_u, fabs(vel->vector->x.data[i][j][k]));
                max_v = fmax(max_v, fabs(vel->vector->y.data[i][j][k]));
                max_w = fmax(max_w, fabs(vel->vector->z.data[i][j][k]));
                max_curvature = fmax(max_curvature, fabs(vel->dummy->curvature[i][j][k]));
            }
        }
    }

    DataType CFL_convection, CFL_gravity, CFL_surface_tension, CFL_diffusion, CFL;
    CFL_convection = CFL_gravity = CFL_surface_tension = CFL_diffusion = CFL = 0.0;

    CFL_convection = max_u / vel->geo->dx + max_v / vel->geo->dy;
    if (vel->vector->x.ndim == 3) {
        CFL_convection += max_w / vel->geo->dz;
    }
    CFL_diffusion = 2.0 / (vel->geo->dx * vel->geo->dx) + 2.0 / (vel->geo->dy * vel->geo->dy);
    if (vel->vector->x.ndim == 3) {
        CFL_diffusion += 2.0 / (vel->geo->dz * vel->geo->dz);
    }
    CFL_diffusion *= vel->params->viscosity_ratio / vel->params->density_ratio / vel->params->Reynolds_number;
    if ( vel->params->Froude_number > 0.0 ){
        CFL_gravity = sqrt(fabs(1.0 - vel->params->density_ratio) / (2 * vel->params->density_ratio)
                           / vel->geo->h / vel->params->Froude_number);
    }
    if (vel->params->Weber_number > 0.0){
        CFL_surface_tension = sqrt(max_curvature / pow(vel->geo->h, 2)
                / vel->params->Weber_number / vel->params->density_ratio);
    }

    CFL = (CFL_convection + CFL_diffusion + sqrt(pow(CFL_convection + CFL_diffusion, 2)
                                                 + 4 * pow(CFL_gravity, 2)
                                                 + 4 * pow(CFL_surface_tension, 2)));

    vel->params->stable_CFL = 1.0 / CFL / vel->geo->h;
    vel->params->dt = fmin(vel->params->stable_CFL, vel->params->max_CFL) * vel->geo->h;
}

void find_density(wrapper *lsf) {
#pragma omp parallel for default(none) shared(lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx; ++i) {
        for (int j = 0; j < lsf->scalar->Ny; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {
                auto h = Heaviside(lsf->scalar->data[i][j][k], lsf->params->ls_width);
                lsf->dummy->density[i][j][k] = h + (1.0 - h) * lsf->params->density_ratio;
            }
        }
    }
}

void find_viscosity(wrapper *lsf) {
#pragma omp parallel for default(none) shared(lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx; ++i) {
        for (int j = 0; j < lsf->scalar->Ny; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {
                auto h = Heaviside(lsf->scalar->data[i][j][k], lsf->params->ls_width);
                lsf->dummy->viscosity[i][j][k] = h + (1.0 - h) * lsf->params->viscosity_ratio;
            }
        }
    }
}
