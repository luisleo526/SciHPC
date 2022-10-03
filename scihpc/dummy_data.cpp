//
// Created by leo on 10/3/22.
//

#include "dummy_data.h"

dummy_data *dummy_data_alloc(scalar_data* f) {

    auto dummy = new dummy_data;

    dummy->grad = init_array(f->Nx, f->Ny, f->Nz);
    dummy->delta = init_array(f->Nx, f->Ny, f->Nz);
    dummy->heaviside = init_array(f->Nx, f->Ny, f->Nz);
    dummy->sign = init_array(f->Nx, f->Ny, f->Nz);
    dummy->tmp = init_array(f->Nx, f->Ny, f->Nz);
    dummy->u_tmp = init_array(f->Nx, f->Ny, f->Nz);
    dummy->v_tmp = init_array(f->Nx, f->Ny, f->Nz);
    dummy->w_tmp = init_array(f->Nx, f->Ny, f->Nz);
    dummy->a = init_array(f->Nx, f->Ny, f->Nz);
    dummy->b = init_array(f->Nx, f->Ny, f->Nz);
    dummy->c = init_array(f->Nx, f->Ny, f->Nz);

    return dummy;
}
