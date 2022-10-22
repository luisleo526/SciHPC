//
// Created by 溫晧良 on 2022/10/1.
//

#include "flux.h"

void identity_flux(scalar_data *f, vector_data *vel) {

#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->data[i][j][k];
            }
        }
    }
}

void identity_flux(scalar_data *f) {

#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->data[i][j][k];
            }
        }
    }

}

void identity_flux(vector_data *vel) {

#pragma omp parallel for default(none) shared(vel) collapse(3)
    for (int i = 0; i < vel->x.Nx; ++i) {
        for (int j = 0; j < vel->x.Ny; ++j) {
            for (int k = 0; k < vel->x.Nz; ++k) {
                vel->x.flux[i][j][k] = vel->x.data[i][j][k];
                vel->y.flux[i][j][k] = vel->y.data[i][j][k];
                vel->z.flux[i][j][k] = vel->z.data[i][j][k];
            }
        }
    }
}


void identity_old_flux(scalar_data *f, vector_data *vel) {

#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->old[i][j][k];
            }
        }
    }
}

void identity_old_flux(scalar_data *f) {

#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->old[i][j][k];
            }
        }
    }

}

void identity_old_flux(vector_data *vel) {

#pragma omp parallel for default(none) shared(vel) collapse(3)
    for (int i = 0; i < vel->x.Nx; ++i) {
        for (int j = 0; j < vel->x.Ny; ++j) {
            for (int k = 0; k < vel->x.Nz; ++k) {
                vel->x.flux[i][j][k] = vel->x.old[i][j][k];
                vel->y.flux[i][j][k] = vel->y.old[i][j][k];
                vel->z.flux[i][j][k] = vel->z.old[i][j][k];
            }
        }
    }
}

void burgers_flux(scalar_data *f, vector_data *vel) {

#pragma omp parallel for default(none) shared(f, vel) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = f->data[i][j][k] * f->data[i][j][k] * 0.5;
            }
        }
    }
}

void identity_with_extrapolation(scalar_data *f) {
    identity_flux(f);

    // x
#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto index_l = f->index_mapping(1, 1, 1);
            auto index_r = f->index_mapping(f->nx, 1, 1);
            auto ir = index_r.i;
            auto il = index_l.i;
            for (int i = 1; i <= f->ghc; ++i) {
                f->flux[ir + i][j][k] = f->flux[ir][j][k];
                f->flux[il - i][j][k] = f->flux[il][j][k];
            }
        }
    }

    if (f->ndim > 1) {
        // y
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, f->ny, 1);
                auto jr = index_r.j;
                auto jl = index_l.j;
                for (int j = 1; j <= f->ghc; ++j) {
                    f->flux[i][jr + j][k] = f->flux[i][jr][k];
                    f->flux[i][jl - j][k] = f->flux[i][jl][k];
                }
            }
        }
    }

    if (f->ndim > 2) {
        // z
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, 1, f->nz);
                auto kr = index_r.k;
                auto kl = index_l.k;
                for (int k = 1; k <= f->ghc; ++k) {
                    f->flux[i][j][kr + k] = f->flux[i][j][kr];
                    f->flux[i][j][kl - k] = f->flux[i][j][kl];
                }
            }
        }
    }
}

void identity_with_extrapolation(scalar_data *f, vector_data *vel) {
    identity_with_extrapolation(f);
}

void identity_with_extrapolation(vector_data *vel) {

    identity_with_extrapolation(&vel->x);
    if (vel->x.ndim > 1) {
        identity_with_extrapolation(&vel->y);
    }
    if (vel->x.ndim > 2) {
        identity_with_extrapolation(&vel->z);
    }
}

void identity_with_extrapolation_face_x(scalar_data *f) {
    identity_flux(f);

    // x
#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto index_l = f->index_mapping(1, 1, 1);
            auto index_r = f->index_mapping(f->nx, 1, 1);
            auto ir = index_r.i;
            auto il = index_l.i;
            for (int i = 1; i <= f->ghc; ++i) {
                f->flux[ir + i][j][k] = f->flux[ir][j][k];
            }
            for (int i = 2; i <= f->ghc; ++i) {
                f->flux[il - i][j][k] = f->flux[il - 1][j][k];
            }
        }
    }

    if (f->ndim > 1) {
        // y
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, f->ny, 1);
                auto jr = index_r.j;
                auto jl = index_l.j;
                for (int j = 1; j <= f->ghc; ++j) {
                    f->flux[i][jr + j][k] = f->flux[i][jr][k];
                    f->flux[i][jl - j][k] = f->flux[i][jl][k];
                }
            }
        }
    }

    if (f->ndim > 2) {
        // z
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, 1, f->nz);
                auto kr = index_r.k;
                auto kl = index_l.k;
                for (int k = 1; k <= f->ghc; ++k) {
                    f->flux[i][j][kr + k] = f->flux[i][j][kr];
                    f->flux[i][j][kl - k] = f->flux[i][j][kl];
                }
            }
        }
    }
}

void identity_with_extrapolation_face_y(scalar_data *f) {
    identity_flux(f);

    // x
#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto index_l = f->index_mapping(1, 1, 1);
            auto index_r = f->index_mapping(f->nx, 1, 1);
            auto ir = index_r.i;
            auto il = index_l.i;
            for (int i = 1; i <= f->ghc; ++i) {
                f->flux[ir + i][j][k] = f->flux[ir][j][k];
                f->flux[il - i][j][k] = f->flux[il][j][k];
            }
        }
    }

    if (f->ndim > 1) {
        // y
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, f->ny, 1);
                auto jr = index_r.j;
                auto jl = index_l.j;
                for (int j = 1; j <= f->ghc; ++j) {
                    f->flux[i][jr + j][k] = f->flux[i][jr][k];
                }
                for (int j = 2; j <= f->ghc; ++j) {
                    f->flux[i][jl - j][k] = f->flux[i][jl - 1][k];
                }
            }
        }
    }

    if (f->ndim > 2) {
        // z
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, 1, f->nz);
                auto kr = index_r.k;
                auto kl = index_l.k;
                for (int k = 1; k <= f->ghc; ++k) {
                    f->flux[i][j][kr + k] = f->flux[i][j][kr];
                    f->flux[i][j][kl - k] = f->flux[i][j][kl];
                }
            }
        }
    }
}

void identity_with_extrapolation_face_z(scalar_data *f) {

    identity_flux(f);

    // x
#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto index_l = f->index_mapping(1, 1, 1);
            auto index_r = f->index_mapping(f->nx, 1, 1);
            auto ir = index_r.i;
            auto il = index_l.i;
            for (int i = 1; i <= f->ghc; ++i) {
                f->flux[ir + i][j][k] = f->flux[ir][j][k];
                f->flux[il - i][j][k] = f->flux[il][j][k];
            }
        }
    }

    if (f->ndim > 1) {
        // y
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, f->ny, 1);
                auto jr = index_r.j;
                auto jl = index_l.j;
                for (int j = 1; j <= f->ghc; ++j) {
                    f->flux[i][jr + j][k] = f->flux[i][jr][k];
                    f->flux[i][jl - j][k] = f->flux[i][jl][k];
                }
            }
        }
    }

    if (f->ndim > 2) {
        // z
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, 1, f->nz);
                auto kr = index_r.k;
                auto kl = index_l.k;
                for (int k = 1; k <= f->ghc; ++k) {
                    f->flux[i][j][kr + k] = f->flux[i][j][kr];
                }
                for (int k = 2; k <= f->ghc; ++k) {
                    f->flux[i][j][kl - k] = f->flux[i][j][kl - 1];
                }
            }
        }
    }
}

void identity_with_extrapolation_face_x(scalar_data *f, vector_data *vel) {
    identity_with_extrapolation_face_x(f);
}

void identity_with_extrapolation_face_y(scalar_data *f, vector_data *vel) {
    identity_with_extrapolation_face_y(f);
}

void identity_with_extrapolation_face_z(scalar_data *f, vector_data *vel) {
    identity_with_extrapolation_face_z(f);
}

void identity_with_extrapolation_face_x(vector_data *f) {
    identity_with_extrapolation_face_x(&f->x);
    if (f->x.ndim > 1) {
        identity_with_extrapolation_face_x(&f->y);
    }
    if (f->x.ndim > 2) {
        identity_with_extrapolation_face_x(&f->z);
    }
}

void identity_with_extrapolation_face_y(vector_data *f) {
    identity_with_extrapolation_face_y(&f->x);
    if (f->x.ndim > 1) {
        identity_with_extrapolation_face_y(&f->y);
    }
    if (f->x.ndim > 2) {
        identity_with_extrapolation_face_y(&f->z);
    }
}

void identity_with_extrapolation_face_z(vector_data *f) {
    identity_with_extrapolation_face_z(&f->x);
    if (f->x.ndim > 1) {
        identity_with_extrapolation_face_z(&f->y);
    }
    if (f->x.ndim > 2) {
        identity_with_extrapolation_face_z(&f->z);
    }
}


//------------------------------------------------------------------

void identity_old_with_extrapolation(scalar_data *f) {
    identity_old_flux(f);

    // x
#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto index_l = f->index_mapping(1, 1, 1);
            auto index_r = f->index_mapping(f->nx, 1, 1);
            auto ir = index_r.i;
            auto il = index_l.i;
            for (int i = 1; i <= f->ghc; ++i) {
                f->flux[ir + i][j][k] = f->flux[ir][j][k];
                f->flux[il - i][j][k] = f->flux[il][j][k];
            }
        }
    }

    if (f->ndim > 1) {
        // y
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, f->ny, 1);
                auto jr = index_r.j;
                auto jl = index_l.j;
                for (int j = 1; j <= f->ghc; ++j) {
                    f->flux[i][jr + j][k] = f->flux[i][jr][k];
                    f->flux[i][jl - j][k] = f->flux[i][jl][k];
                }
            }
        }
    }

    if (f->ndim > 2) {
        // z
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, 1, f->nz);
                auto kr = index_r.k;
                auto kl = index_l.k;
                for (int k = 1; k <= f->ghc; ++k) {
                    f->flux[i][j][kr + k] = f->flux[i][j][kr];
                    f->flux[i][j][kl - k] = f->flux[i][j][kl];
                }
            }
        }
    }
}

void identity_old_with_extrapolation(scalar_data *f, vector_data *vel) {
    identity_old_with_extrapolation(f);
}

void identity_old_with_extrapolation(vector_data *vel) {

    identity_old_with_extrapolation(&vel->x);
    if (vel->x.ndim > 1) {
        identity_old_with_extrapolation(&vel->y);
    }
    if (vel->x.ndim > 2) {
        identity_old_with_extrapolation(&vel->z);
    }
}

void identity_old_with_extrapolation_face_x(scalar_data *f) {
    identity_old_flux(f);

    // x
#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto index_l = f->index_mapping(1, 1, 1);
            auto index_r = f->index_mapping(f->nx, 1, 1);
            auto ir = index_r.i;
            auto il = index_l.i;
            for (int i = 1; i <= f->ghc; ++i) {
                f->flux[ir + i][j][k] = f->flux[ir][j][k];
            }
            for (int i = 2; i <= f->ghc; ++i) {
                f->flux[il - i][j][k] = f->flux[il - 1][j][k];
            }
        }
    }

    if (f->ndim > 1) {
        // y
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, f->ny, 1);
                auto jr = index_r.j;
                auto jl = index_l.j;
                for (int j = 1; j <= f->ghc; ++j) {
                    f->flux[i][jr + j][k] = f->flux[i][jr][k];
                    f->flux[i][jl - j][k] = f->flux[i][jl][k];
                }
            }
        }
    }

    if (f->ndim > 2) {
        // z
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, 1, f->nz);
                auto kr = index_r.k;
                auto kl = index_l.k;
                for (int k = 1; k <= f->ghc; ++k) {
                    f->flux[i][j][kr + k] = f->flux[i][j][kr];
                    f->flux[i][j][kl - k] = f->flux[i][j][kl];
                }
            }
        }
    }
}

void identity_old_with_extrapolation_face_y(scalar_data *f) {
    identity_old_flux(f);

    // x
#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto index_l = f->index_mapping(1, 1, 1);
            auto index_r = f->index_mapping(f->nx, 1, 1);
            auto ir = index_r.i;
            auto il = index_l.i;
            for (int i = 1; i <= f->ghc; ++i) {
                f->flux[ir + i][j][k] = f->flux[ir][j][k];
                f->flux[il - i][j][k] = f->flux[il][j][k];
            }
        }
    }

    if (f->ndim > 1) {
        // y
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, f->ny, 1);
                auto jr = index_r.j;
                auto jl = index_l.j;
                for (int j = 1; j <= f->ghc; ++j) {
                    f->flux[i][jr + j][k] = f->flux[i][jr][k];
                }
                for (int j = 2; j <= f->ghc; ++j) {
                    f->flux[i][jl - j][k] = f->flux[i][jl - 1][k];
                }
            }
        }
    }

    if (f->ndim > 2) {
        // z
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, 1, f->nz);
                auto kr = index_r.k;
                auto kl = index_l.k;
                for (int k = 1; k <= f->ghc; ++k) {
                    f->flux[i][j][kr + k] = f->flux[i][j][kr];
                    f->flux[i][j][kl - k] = f->flux[i][j][kl];
                }
            }
        }
    }
}

void identity_old_with_extrapolation_face_z(scalar_data *f) {

    identity_old_flux(f);

    // x
#pragma omp parallel for default(none) shared(f) collapse(2)
    for (int j = 0; j < f->Ny; ++j) {
        for (int k = 0; k < f->Nz; ++k) {
            auto index_l = f->index_mapping(1, 1, 1);
            auto index_r = f->index_mapping(f->nx, 1, 1);
            auto ir = index_r.i;
            auto il = index_l.i;
            for (int i = 1; i <= f->ghc; ++i) {
                f->flux[ir + i][j][k] = f->flux[ir][j][k];
                f->flux[il - i][j][k] = f->flux[il][j][k];
            }
        }
    }

    if (f->ndim > 1) {
        // y
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int k = 0; k < f->Nz; ++k) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, f->ny, 1);
                auto jr = index_r.j;
                auto jl = index_l.j;
                for (int j = 1; j <= f->ghc; ++j) {
                    f->flux[i][jr + j][k] = f->flux[i][jr][k];
                    f->flux[i][jl - j][k] = f->flux[i][jl][k];
                }
            }
        }
    }

    if (f->ndim > 2) {
        // z
#pragma omp parallel for default(none) shared(f) collapse(2)
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                auto index_l = f->index_mapping(1, 1, 1);
                auto index_r = f->index_mapping(1, 1, f->nz);
                auto kr = index_r.k;
                auto kl = index_l.k;
                for (int k = 1; k <= f->ghc; ++k) {
                    f->flux[i][j][kr + k] = f->flux[i][j][kr];
                }
                for (int k = 2; k <= f->ghc; ++k) {
                    f->flux[i][j][kl - k] = f->flux[i][j][kl - 1];
                }
            }
        }
    }
}

void identity_old_with_extrapolation_face_x(scalar_data *f, vector_data *vel) {
    identity_old_with_extrapolation_face_x(f);
}

void identity_old_with_extrapolation_face_y(scalar_data *f, vector_data *vel) {
    identity_old_with_extrapolation_face_y(f);
}

void identity_old_with_extrapolation_face_z(scalar_data *f, vector_data *vel) {
    identity_old_with_extrapolation_face_z(f);
}

void identity_old_with_extrapolation_face_x(vector_data *f) {
    identity_old_with_extrapolation_face_x(&f->x);
    if (f->x.ndim > 1) {
        identity_old_with_extrapolation_face_x(&f->y);
    }
    if (f->x.ndim > 2) {
        identity_old_with_extrapolation_face_x(&f->z);
    }
}

void identity_old_with_extrapolation_face_y(vector_data *f) {
    identity_old_with_extrapolation_face_y(&f->x);
    if (f->x.ndim > 1) {
        identity_old_with_extrapolation_face_y(&f->y);
    }
    if (f->x.ndim > 2) {
        identity_old_with_extrapolation_face_y(&f->z);
    }
}

void identity_old_with_extrapolation_face_z(vector_data *f) {
    identity_old_with_extrapolation_face_z(&f->x);
    if (f->x.ndim > 1) {
        identity_old_with_extrapolation_face_z(&f->y);
    }
    if (f->x.ndim > 2) {
        identity_old_with_extrapolation_face_z(&f->z);
    }
}

void guess_from_history(scalar_data *f) {
#pragma omp parallel for default(none) shared(f) collapse(3)
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                f->flux[i][j][k] = 2.0 * f->old[i][j][k] - f->old2[i][j][k];
            }
        }
    }
}
