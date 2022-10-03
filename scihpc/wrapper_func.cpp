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
                if (lsf->scalar->ndim > 2){
                    gradient += pow(lsf->scalar->fz[i][j][k], 2);
                }
                lsf->dummy->grad[i][j][k] = sqrt(gradient);
            }
        }
    }
}

void store(wrapper *f) {
    if (f->is_scalar){
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    f->dummy->tmp[i][j][k] = f->scalar->data[i][j][k];
                }
            }
        }
    }else{
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
