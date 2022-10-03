//
// Created by leo on 10/3/22.
//

#include "dummy_data.h"

dummy_data *dummy_data_alloc(int nx, int ny, int nz) {
    auto dummy = new dummy_data;
    dummy->heaviside = init_array(nx, ny, nz);
    dummy->sign = init_array(nx, ny, nz);
    dummy->delta = init_array(nx, ny, nz);
    dummy->tmp = init_array(nx, ny, nz);
    dummy->u_tmp = init_array(nx, ny, nz);
    dummy->v_tmp = init_array(nx, ny, nz);
    dummy->w_tmp = init_array(nx, ny, nz);
    dummy->a = init_array(nx, ny, nz);
    dummy->b = init_array(nx, ny, nz);
    dummy->c = init_array(nx, ny, nz);
    return dummy;
}
