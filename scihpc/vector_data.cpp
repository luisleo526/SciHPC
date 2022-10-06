//
// Created by 溫晧良 on 2022/9/28.
//

#include "vector_data.h"

vector_data::vector_data(const int _nx) :
        x(_nx), y(_nx), z(_nx) {
    x = scalar_data(_nx);
    y = scalar_data(_nx);
    z = scalar_data(_nx);
}

vector_data::vector_data(const int _nx, const int _ny) :
        x(_nx, _ny), y(_nx, _ny), z(_nx, _ny) {
    x = scalar_data(_nx, _ny);
    y = scalar_data(_nx, _ny);
    z = scalar_data(_nx, _ny);
}

vector_data::vector_data(const int _nx, const int _ny, const int _nz) :
        x(_nx, _ny, _nz), y(_nx, _ny, _nz), z(_nx, _ny, _nz) {
    x = scalar_data(_nx, _ny, _nz);
    y = scalar_data(_nx, _ny, _nz);
    z = scalar_data(_nx, _ny, _nz);
}

void vector_data::store() {
    x.store();
    y.store();
    z.store();
}
