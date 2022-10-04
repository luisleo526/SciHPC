//
// Created by 溫晧良 on 2022/10/5.
//

#include "projection_method.h"

projection_method::projection_method(scalar_data *f) {
    u_src= init_array(f->Nx, f->Ny, f->Nz);
    v_src= init_array(f->Nx, f->Ny, f->Nz);
    w_src= init_array(f->Nx, f->Ny, f->Nz);

    u_src_old= init_array(f->Nx, f->Ny, f->Nz);
    v_src_old= init_array(f->Nx, f->Ny, f->Nz);
    w_src_old= init_array(f->Nx, f->Ny, f->Nz);
}
