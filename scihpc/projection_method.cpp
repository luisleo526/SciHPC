//
// Created by 溫晧良 on 2022/10/5.
//

#include "projection_method.h"
#include <iostream>

projection_method::projection_method(scalar_data *f) {
    u_src = init_array(f->Nx, f->Ny, f->Nz);
    v_src = init_array(f->Nx, f->Ny, f->Nz);
    w_src = init_array(f->Nx, f->Ny, f->Nz);

    u_src_old = init_array(f->Nx, f->Ny, f->Nz);
    v_src_old = init_array(f->Nx, f->Ny, f->Nz);
    w_src_old = init_array(f->Nx, f->Ny, f->Nz);

    // For PPE iterative solver
    CC = init_array(f->Nx, f->Ny, f->Nz);
    CR = init_array(f->Nx, f->Ny, f->Nz);
    CL = init_array(f->Nx, f->Ny, f->Nz);
    CU = init_array(f->Nx, f->Ny, f->Nz);
    CD = init_array(f->Nx, f->Ny, f->Nz);
    CF = init_array(f->Nx, f->Ny, f->Nz);
    CB = init_array(f->Nx, f->Ny, f->Nz);
    RHS = init_array(f->Nx, f->Ny, f->Nz);
}

void projection_method::add_stress_x(wrapper *nvel, wrapper *lsf) {

#pragma omp parallel for default(none) shared(nvel, lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx - 1; ++i) {
        for (int j = 0; j < lsf->scalar->Ny; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {

                auto h = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i + 1][j][k]);
                auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i + 1][j][k]);
                auto mu = h + (1.0 - h) * lsf->params->viscosity_ratio;
                auto rho = h + (1.0 - h) * lsf->params->density_ratio;

                auto phix = 0.5 * (lsf->scalar->fx[i][j][k] + lsf->scalar->fx[i + 1][j][k]);
                auto phiy = 0.5 * (lsf->scalar->fy[i][j][k] + lsf->scalar->fy[i + 1][j][k]);
                auto phiz = 0.5 * (lsf->scalar->fz[i][j][k] + lsf->scalar->fz[i + 1][j][k]);

                auto ux = nvel->vector->x.fx[i][j][k];
                auto uy = nvel->vector->x.fy[i][j][k];
                auto uz = nvel->vector->x.fz[i][j][k];

                auto uxx = nvel->vector->x.fxx[i][j][k];
                auto uyy = nvel->vector->x.fyy[i][j][k];
                auto uzz = nvel->vector->x.fzz[i][j][k];

                auto vx = nvel->vector->y.fx[i][j][k];
                auto wx = nvel->vector->z.fx[i][j][k];

                auto stress_part1 = uxx + uyy;
                auto stress_part2 = 2.0 * phix * ux + phiy * (uy + vx);
                if (lsf->scalar->ndim > 2) {
                    stress_part1 += uzz;
                    stress_part2 += phiz * (uz + wx);
                }
                stress_part1 *= mu;
                stress_part2 *= delta * (1.0 - lsf->params->viscosity_ratio);

                u_src[i][j][k] += (stress_part1 + stress_part2) / lsf->params->Reynolds_number / rho;
            }
        }
    }

}

void projection_method::add_stress_y(wrapper *nvel, wrapper *lsf) {

#pragma omp parallel for default(none) shared(nvel, lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx; ++i) {
        for (int j = 0; j < lsf->scalar->Ny - 1; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {

                auto h = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j + 1][k]);
                auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i][j + 1][k]);
                auto mu = h + (1.0 - h) * lsf->params->viscosity_ratio;
                auto rho = h + (1.0 - h) * lsf->params->density_ratio;

                auto phix = 0.5 * (lsf->scalar->fx[i][j][k] + lsf->scalar->fx[i][j + 1][k]);
                auto phiy = 0.5 * (lsf->scalar->fy[i][j][k] + lsf->scalar->fy[i][j + 1][k]);
                auto phiz = 0.5 * (lsf->scalar->fz[i][j][k] + lsf->scalar->fz[i][j + 1][k]);

                auto vx = nvel->vector->y.fx[i][j][k];
                auto vy = nvel->vector->y.fy[i][j][k];
                auto vz = nvel->vector->y.fz[i][j][k];

                auto vxx = nvel->vector->y.fxx[i][j][k];
                auto vyy = nvel->vector->y.fyy[i][j][k];
                auto vzz = nvel->vector->y.fzz[i][j][k];

                auto uy = nvel->vector->x.fy[i][j][k];
                auto wy = nvel->vector->z.fy[i][j][k];

                auto stress_part1 = vxx + vyy;
                auto stress_part2 = 2.0 * phiy * vy + phix * (vx + uy);
                if (lsf->scalar->ndim > 2) {
                    stress_part1 += vzz;
                    stress_part2 += phiz * (wy + vz);
                }
                stress_part1 *= mu;
                stress_part2 *= delta * (1.0 - lsf->params->viscosity_ratio);

                v_src[i][j][k] += (stress_part1 + stress_part2) / lsf->params->Reynolds_number / rho;
            }
        }
    }
}

void projection_method::add_stress_z(wrapper *nvel, wrapper *lsf) {

#pragma omp parallel for default(none) shared(nvel, lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx; ++i) {
        for (int j = 0; j < lsf->scalar->Ny; ++j) {
            for (int k = 0; k < lsf->scalar->Nz - 1; ++k) {

                auto h = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j][k + 1]);
                auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i][j][k + 1]);
                auto mu = h + (1.0 - h) * lsf->params->viscosity_ratio;
                auto rho = h + (1.0 - h) * lsf->params->density_ratio;

                auto phix = 0.5 * (lsf->scalar->fx[i][j][k] + lsf->scalar->fx[i][j][k + 1]);
                auto phiy = 0.5 * (lsf->scalar->fy[i][j][k] + lsf->scalar->fy[i][j][k + 1]);
                auto phiz = 0.5 * (lsf->scalar->fz[i][j][k] + lsf->scalar->fz[i][j][k + 1]);

                auto wx = nvel->vector->z.fx[i][j][k];
                auto wy = nvel->vector->z.fy[i][j][k];
                auto wz = nvel->vector->z.fz[i][j][k];

                auto wxx = nvel->vector->z.fxx[i][j][k];
                auto wyy = nvel->vector->z.fyy[i][j][k];
                auto wzz = nvel->vector->z.fzz[i][j][k];

                auto uz = nvel->vector->x.fz[i][j][k];
                auto vz = nvel->vector->y.fz[i][j][k];

                auto stress_part1 = wxx + wyy + wzz;
                auto stress_part2 = 2.0 * phiz * wz + phix * (wx + uz) + phiy * (wy + vz);

                stress_part1 *= mu;
                stress_part2 *= delta * (1.0 - lsf->params->viscosity_ratio);

                w_src[i][j][k] += (stress_part1 + stress_part2) / lsf->params->Reynolds_number / rho;
            }
        }
    }
}

void projection_method::find_source(wrapper *vel, wrapper *nvel, wrapper *lsf, structured_grid *geo) {

#pragma omp parallel for default(none) collapse(3)
    for (int i = 0; i < vel->vector->x.Nx; ++i) {
        for (int j = 0; j < vel->vector->x.Ny; ++j) {
            for (int k = 0; k < vel->vector->x.Nz; ++k) {
                u_src_old[i][j][k] = u_src[i][j][k];
                v_src_old[i][j][k] = v_src[i][j][k];
                w_src_old[i][j][k] = w_src[i][j][k];
            }
        }
    }

    // prepare level set function
    find_delta(lsf);
    find_heavyside(lsf);
    identity_flux(lsf->scalar);
    find_curvature(lsf);

    // x direction
    all_to_face_x(vel, nvel);
    no_slip_face_x(nvel->vector);
    identity_flux(nvel->vector);
    vel->scalar = &vel->vector->x;
    convection(vel, nvel, geo, u_src, &identity_flux);
    nvel->solvers->ccd->find_derivatives(nvel->vector);
    add_stress_x(nvel, lsf);

    // y direction
    all_to_face_y(vel, nvel);
    no_slip_face_y(nvel->vector);
    identity_flux(nvel->vector);
    vel->scalar = &vel->vector->y;
    convection(vel, nvel, geo, v_src, &identity_flux);
    nvel->solvers->ccd->find_derivatives(nvel->vector);
    add_stress_y(nvel, lsf);

    // z direction
    if (lsf->scalar->ndim > 2) {
        all_to_face_z(vel, nvel);
        no_slip_face_z(nvel->vector);
        identity_flux(nvel->vector);
        vel->scalar = &vel->vector->z;
        convection(vel, nvel, geo, w_src, &identity_flux);
        nvel->solvers->ccd->find_derivatives(nvel->vector);
        add_stress_z(nvel, lsf);
    }

    // add gravity
    if (lsf->params->Froude_number > 0.0) {
        if (lsf->scalar->ndim == 2) {
#pragma omp parallel for default(none) shared(nvel, lsf) collapse(3)
            for (int i = 0; i < lsf->scalar->Nx; ++i) {
                for (int j = 0; j < lsf->scalar->Ny; ++j) {
                    for (int k = 0; k < lsf->scalar->Nz; ++k) {
                        v_src[i][j][k] += 1.0 / lsf->params->Froude_number;
                    }
                }
            }
        } else {
#pragma omp parallel for default(none) shared(nvel, lsf) collapse(3)
            for (int i = 0; i < lsf->scalar->Nx; ++i) {
                for (int j = 0; j < lsf->scalar->Ny; ++j) {
                    for (int k = 0; k < lsf->scalar->Nz; ++k) {
                        w_src[i][j][k] += 1.0 / lsf->params->Froude_number;
                    }
                }
            }
        }
    }

    if (lsf->params->Weber_number > 0.0) {
#pragma omp parallel for default(none) shared(nvel, lsf) collapse(3)
        for (int i = 0; i < lsf->scalar->Nx - 1; ++i) {
            for (int j = 0; j < lsf->scalar->Ny - 1; ++j) {
                for (int k = 0; k < lsf->scalar->Nz; ++k) {
                    // x direction
                    auto h = 0.5 * (lsf->dummy->heaviside[i + 1][j][k] + lsf->dummy->heaviside[i][j][k]);
                    auto delta = 0.5 * (lsf->dummy->delta[i + 1][j][k] + lsf->dummy->delta[i][j][k]);
                    auto curv = 0.5 * (lsf->dummy->curvature[i + 1][j][k] + lsf->dummy->curvature[i][j][k]);
                    auto dir = 0.5 * (lsf->scalar->fx[i + 1][j][k] + lsf->scalar->fx[i][j][k]);
                    auto density = h + (1.0 - h) * lsf->params->density_ratio;
                    u_src[i][j][k] -= curv * dir * delta / (density * lsf->params->Weber_number);

                    // y direction
                    h = 0.5 * (lsf->dummy->heaviside[i][j + 1][k] + lsf->dummy->heaviside[i][j][k]);
                    delta = 0.5 * (lsf->dummy->delta[i][j + 1][k] + lsf->dummy->delta[i][j][k]);
                    curv = 0.5 * (lsf->dummy->curvature[i][j + 1][k] + lsf->dummy->curvature[i][j][k]);
                    dir = 0.5 * (lsf->scalar->fy[i][j + 1][k] + lsf->scalar->fy[i][j][k]);
                    density = h + (1.0 - h) * lsf->params->density_ratio;
                    v_src[i][j][k] -= curv * dir * delta / (density * lsf->params->Weber_number);
                }
            }
        }

        if (lsf->scalar->ndim > 2) {
#pragma omp parallel for default(none) shared(nvel, lsf) collapse(3)
            for (int i = 0; i < lsf->scalar->Nx; ++i) {
                for (int j = 0; j < lsf->scalar->Ny; ++j) {
                    for (int k = 0; k < lsf->scalar->Nz - 1; ++k) {
                        // z direction
                        auto h = 0.5 * (lsf->dummy->heaviside[i][j][k + 1] + lsf->dummy->heaviside[i][j][k]);
                        auto delta = 0.5 * (lsf->dummy->delta[i][j][k + 1] + lsf->dummy->delta[i][j][k]);
                        auto curv = 0.5 * (lsf->dummy->curvature[i][j][k + 1] + lsf->dummy->curvature[i][j][k]);
                        auto dir = 0.5 * (lsf->scalar->fz[i][j][k + 1] + lsf->scalar->fz[i][j][k]);
                        auto density = h + (1.0 - h) * lsf->params->density_ratio;
                        w_src[i][j][k] -= curv * dir * delta / (density * lsf->params->Weber_number);
                    }
                }
            }
        }

    }

    vel->scalar = nullptr;
}

void projection_method::find_intermediate_velocity(wrapper *vel) {

#pragma omp parallel for default(none) shared(vel) collapse(3)
    for (int i = 0; i < vel->vector->x.Nx; ++i) {
        for (int j = 0; j < vel->vector->x.Ny; ++j) {
            for (int k = 0; k < vel->vector->x.Nz; ++k) {
                vel->vector->x.data[i][j][k] += vel->params->dt * (1.5 * u_src[i][j][k] - 0.5 * u_src_old[i][j][k]);
                vel->vector->y.data[i][j][k] += vel->params->dt * (1.5 * v_src[i][j][k] - 0.5 * v_src_old[i][j][k]);
                vel->vector->z.data[i][j][k] += vel->params->dt * (1.5 * w_src[i][j][k] - 0.5 * w_src_old[i][j][k]);
            }
        }
    }

    no_slip_face(vel->vector);
}

void projection_method::solve_ppe(wrapper *pressure, wrapper *lsf, wrapper *vel, structured_grid *geo) {

    if (lsf->scalar->ndim == 2) {
#pragma omp parallel for default(none) shared(lsf, geo) collapse(3)
        for (int i = 1; i < lsf->scalar->Nx - 1; ++i) {
            for (int j = 1; j < lsf->scalar->Ny - 1; ++j) {
                for (int k = 0; k < lsf->scalar->Nz; ++k) {
                    auto R = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i + 1][j][k]);
                    auto L = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i - 1][j][k]);
                    auto U = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j + 1][k]);
                    auto D = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j - 1][k]);

                    R = R * (1.0 - R) * lsf->params->density_ratio;
                    L = L * (1.0 - L) * lsf->params->density_ratio;
                    U = U * (1.0 - U) * lsf->params->density_ratio;
                    D = D * (1.0 - D) * lsf->params->density_ratio;

                    CR[i][j][k] = 1.0 / R / geo->dx / geo->dx;
                    CL[i][j][k] = 1.0 / L / geo->dx / geo->dx;
                    CU[i][j][k] = 1.0 / U / geo->dy / geo->dy;
                    CD[i][j][k] = 1.0 / D / geo->dy / geo->dy;

                    if (i == lsf->scalar->ghc) {
                        CL[i][j][k] = 0.0;
                    }

                    if (i == lsf->scalar->Nx - lsf->scalar->ghc - 1) {
                        CR[i][j][k] = 0.0;
                    }

                    if (j == lsf->scalar->ghc) {
                        CD[i][j][k] = 0.0;
                    }

                    if (j == lsf->scalar->Ny - lsf->scalar->ghc - 1) {
                        CU[i][j][k] = 0.0;
                    }

                    CC[i][j][k] = -(CR[i][j][k] + CL[i][j][k] + CU[i][j][k] + CD[i][j][k]);

                    RHS[i][j][k] = (vel->vector->x.data[i][j][k] - vel->vector->x.data[i - 1][j][k]) / geo->dx +
                                   (vel->vector->y.data[i][j][k] - vel->vector->y.data[i][j - 1][k]) / geo->dy;
                    RHS[i][j][k] /= lsf->params->dt;
                }
            }
        }
    } else if (lsf->scalar->ndim == 3) {
#pragma omp parallel for default(none) shared(lsf, geo) collapse(3)
        for (int i = 1; i < lsf->scalar->Nx - 1; ++i) {
            for (int j = 1; j < lsf->scalar->Ny - 1; ++j) {
                for (int k = 1; k < lsf->scalar->Nz - 1; ++k) {
                    auto R = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i + 1][j][k]);
                    auto L = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i - 1][j][k]);
                    auto U = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j + 1][k]);
                    auto D = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j - 1][k]);
                    auto F = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j][k + 1]);
                    auto B = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j][k - 1]);

                    R = R * (1.0 - R) * lsf->params->density_ratio;
                    L = L * (1.0 - L) * lsf->params->density_ratio;
                    U = U * (1.0 - U) * lsf->params->density_ratio;
                    D = D * (1.0 - D) * lsf->params->density_ratio;
                    F = F * (1.0 - F) * lsf->params->density_ratio;
                    B = B * (1.0 - B) * lsf->params->density_ratio;

                    CR[i][j][k] = 1.0 / R / geo->dx / geo->dx;
                    CL[i][j][k] = 1.0 / L / geo->dx / geo->dx;
                    CU[i][j][k] = 1.0 / U / geo->dy / geo->dy;
                    CD[i][j][k] = 1.0 / D / geo->dy / geo->dy;
                    CF[i][j][k] = 1.0 / F / geo->dz / geo->dz;
                    CB[i][j][k] = 1.0 / B / geo->dz / geo->dz;

                    if (i == lsf->scalar->ghc) {
                        CL[i][j][k] = 0.0;
                    }

                    if (i == lsf->scalar->Nx - lsf->scalar->ghc - 1) {
                        CR[i][j][k] = 0.0;
                    }

                    if (j == lsf->scalar->ghc) {
                        CD[i][j][k] = 0.0;
                    }

                    if (j == lsf->scalar->Ny - lsf->scalar->ghc - 1) {
                        CU[i][j][k] = 0.0;
                    }

                    if (k == lsf->scalar->ghc) {
                        CB[i][j][k] = 0.0;
                    }

                    if (k == lsf->scalar->Nz - lsf->scalar->ghc - 1) {
                        CF[i][j][k] = 0.0;
                    }

                    CC[i][j][k] = -(CR[i][j][k] + CL[i][j][k] + CU[i][j][k] + CD[i][j][k] + CF[i][j][k] + CB[i][j][k]);

                    RHS[i][j][k] = (vel->vector->x.data[i][j][k] - vel->vector->x.data[i - 1][j][k]) / geo->dx +
                                   (vel->vector->y.data[i][j][k] - vel->vector->y.data[i][j - 1][k]) / geo->dy +
                                   (vel->vector->z.data[i][j][k] - vel->vector->z.data[i][j][k - 1]) / geo->dz;
                    RHS[i][j][k] /= vel->params->dt;
                }
            }
        }
    }

    // solve pressure poisson equation

    int iter = 0;
    DataType error, sump;
    do {
        store_tmp(pressure);
        if (pressure->scalar->ndim == 2) {
            sump = 0.0;
#pragma omp parallel for default(none) shared(lsf, geo) collapse(2) reduction(+:sump)
            for (int i = 0; i < pressure->scalar->nx; ++i) {
                for (int j = 0; j < pressure->scalar->ny; ++j) {
                    auto index = pressure->scalar->index_mapping(i + 1, j + 1, 1);
                    auto I = index.i;
                    auto J = index.j;
                    auto K = index.k;

                    auto pnew = RHS[I][J][K]
                                - (CL[I][J][K] * pressure->scalar->data[I - 1][J][K] +
                                   CR[I][J][K] * pressure->scalar->data[I + 1][J][K] +
                                   CD[I][J][K] * pressure->scalar->data[I][J - 1][K] +
                                   CU[I][J][K] * pressure->scalar->data[I][J + 1][K]) / CC[I][J][K];
                    pressure->scalar->data[I][J][K] = 1.5 * pnew - 0.5 * pressure->scalar->data[I][J][K];
                    sump += pressure->scalar->data[I][J][K];
                }
            }
            sump = sump / (pressure->scalar->nx * pressure->scalar->ny);
#pragma omp parallel for default(none) shared(lsf, geo) collapse(2)
            for (int i = 0; i < pressure->scalar->nx; ++i) {
                for (int j = 0; j < pressure->scalar->ny; ++j) {
                    auto index = pressure->scalar->index_mapping(i + 1, j + 1, 1);
                    pressure->scalar->data[index.i][index.j][index.k] -= sump;
                }
            }

        } else {
#pragma omp parallel for default(none) shared(lsf, geo) collapse(3)
            for (int i = 0; i < pressure->scalar->nx; ++i) {
                for (int j = 0; j < pressure->scalar->ny; ++j) {
                    for (int k = 0; k < pressure->scalar->nz; ++k) {
                        auto index = pressure->scalar->index_mapping(i + 1, j + 1, k + 1);
                        auto I = index.i;
                        auto J = index.j;
                        auto K = index.k;
                        auto pnew = RHS[I][J][K]
                                    - (CL[I][J][K] * pressure->scalar->data[I - 1][J][K] +
                                       CR[I][J][K] * pressure->scalar->data[I + 1][J][K] +
                                       CD[I][J][K] * pressure->scalar->data[I][J - 1][K] +
                                       CU[I][J][K] * pressure->scalar->data[I][J + 1][K] +
                                       CB[I][J][K] * pressure->scalar->data[I][J][K - 1] +
                                       CF[I][J][K] * pressure->scalar->data[I][J][K + 1]) / CC[I][J][K];

                        pressure->scalar->data[I][J][K] = 1.5 * pnew - 0.5 * pressure->scalar->data[I][J][K];
                        sump += pressure->scalar->data[I][J][K];
                    }
                }
            }
        }
        zero_order_extrapolation(pressure->scalar);
        error = l2norm(pressure);
        if (iter % 1000 == 0) {
            std::cout << "iter: " << iter << " error: " << error << std::endl;
        }
    } while (++iter < pressure->params->ppe_max_iter and error > pressure->params->ppe_tol);

    std::cout << "iter: " << iter << " error: " << error << std::endl;

}

void projection_method::find_final_velocity(wrapper *vel, wrapper *pressure, wrapper *lsf, structured_grid *geo) {

    if (vel->vector->x.ndim == 2) {
#pragma omp parallel for default(none) shared(vel, pressure, geo) collapse(2)
        for (int i = 0; i < vel->vector->x.Nx - 1; ++i) {
            for (int j = 0; j < vel->vector->x.Ny - 1; ++j) {
                auto h = 0.5 * (lsf->dummy->heaviside[i][j][0] + lsf->dummy->heaviside[i + 1][j][0]);
                auto density = (1 - h) * vel->params->density_ratio + h;
                vel->vector->x.data[i][j][0] -= pressure->params->dt *
                                                (pressure->scalar->data[i + 1][j][0] -
                                                 pressure->scalar->data[i][j][0]) / geo->dx / density;

                h = 0.5 * (lsf->dummy->heaviside[i][j][0] + lsf->dummy->heaviside[i][j + 1][0]);
                density = (1 - h) * vel->params->density_ratio + h;

                vel->vector->y.data[i][j][0] -= pressure->params->dt *
                                                (pressure->scalar->data[i][j + 1][0] -
                                                 pressure->scalar->data[i][j][0]) / geo->dy / density;
            }
        }
    } else {
#pragma omp parallel for default(none) shared(vel, pressure, geo) collapse(3)
        for (int i = 0; i < vel->vector->x.Nx - 1; ++i) {
            for (int j = 0; j < vel->vector->x.Ny - 1; ++j) {
                for (int k = 0; k < vel->vector->x.Nz - 1; ++k) {
                    auto h = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i + 1][j][k]);
                    auto density = (1 - h) * vel->params->density_ratio + h;
                    vel->vector->x.data[i][j][k] -= pressure->params->dt *
                                                    (pressure->scalar->data[i + 1][j][k] -
                                                     pressure->scalar->data[i][j][k]) / geo->dx / density;

                    h = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j + 1][k]);
                    density = (1 - h) * vel->params->density_ratio + h;
                    vel->vector->y.data[i][j][k] -= pressure->params->dt *
                                                    (pressure->scalar->data[i][j + 1][k] -
                                                     pressure->scalar->data[i][j][k]) / geo->dy / density;

                    h = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j][k + 1]);
                    density = (1 - h) * vel->params->density_ratio + h;
                    vel->vector->z.data[i][j][k] -= pressure->params->dt *
                                                    (pressure->scalar->data[i][j][k + 1] -
                                                     pressure->scalar->data[i][j][k]) / geo->dz / density;
                }
            }
        }
    }

    no_slip_face(vel->vector);
}

void projection_method::solve(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf, structured_grid *geo) {
    find_source(vel, nvel, lsf, geo);
    find_intermediate_velocity(vel);
    solve_ppe(pressure, lsf, vel, geo);
    find_final_velocity(vel, pressure, lsf, geo);
    node_from_face(vel, nvel);
    no_slip_node(nvel->vector);
}
