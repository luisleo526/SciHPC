//
// Created by 溫晧良 on 2022/9/28.
//

#include "scalar_data.h"
#include <cassert>

scalar_data::scalar_data(const int _nx) {
    ndim = 1;
    nx = _nx;
    ny = 1;
    nz = 1;

    Nx = _nx + 2 * ghc;
    Ny = 1;
    Nz = 1;

    data = init_array(Nx, Ny, Nz);

    fx = init_array(Nx, Ny, Nz);
    fy = init_array(Nx, Ny, Nz);
    fz = init_array(Nx, Ny, Nz);

    fxx = init_array(Nx, Ny, Nz);
    fyy = init_array(Nx, Ny, Nz);
    fzz = init_array(Nx, Ny, Nz);

    fxy = init_array(Nx, Ny, Nz);
    fyz = init_array(Nx, Ny, Nz);
    fxz = init_array(Nx, Ny, Nz);

}

scalar_data::scalar_data(const int _nx, const int _ny) {
    ndim = 2;
    nx = _nx;
    ny = _ny;
    nz = 1;

    Nx = _nx + 2 * ghc;
    Ny = _ny + 2 * ghc;
    Nz = 1;

    data = init_array(Nx, Ny, Nz);

    fx = init_array(Nx, Ny, Nz);
    fy = init_array(Nx, Ny, Nz);
    fz = init_array(Nx, Ny, Nz);

    fxx = init_array(Nx, Ny, Nz);
    fyy = init_array(Nx, Ny, Nz);
    fzz = init_array(Nx, Ny, Nz);

    fxy = init_array(Nx, Ny, Nz);
    fyz = init_array(Nx, Ny, Nz);
    fxz = init_array(Nx, Ny, Nz);

}

scalar_data::scalar_data(const int _nx, const int _ny, const int _nz) {
    ndim = 3;
    nx = _nx;
    ny = _ny;
    nz = _nz;

    Nx = _nx + 2 * ghc;
    Ny = _ny + 2 * ghc;
    Nz = _nz + 2 * ghc;

    data = init_array(Nx, Ny, Nz);

    fx = init_array(Nx, Ny, Nz);
    fy = init_array(Nx, Ny, Nz);
    fz = init_array(Nx, Ny, Nz);

    fxx = init_array(Nx, Ny, Nz);
    fyy = init_array(Nx, Ny, Nz);
    fzz = init_array(Nx, Ny, Nz);

    fxy = init_array(Nx, Ny, Nz);
    fyz = init_array(Nx, Ny, Nz);
    fxz = init_array(Nx, Ny, Nz);
}

indices scalar_data::index_mapping(int i, int j, int k) {
    int I, J, K;

    bool i_checked = (nx > 1 and i > -ghc and i < nx + ghc + 1) or nx == i;
    bool j_checked = (ny > 1 and j > -ghc and j < ny + ghc + 1) or ny == j;
    bool k_checked = (nz > 1 and k > -ghc and k < nz + ghc + 1) or nz == k;

    assert(i_checked);
    assert(j_checked);
    assert(k_checked);

    I = i - 1;
    J = j - 1;
    K = k - 1;

    if (nx > 1) {
        I = I + ghc;
    }

    if (ny > 1) {
        J = J + ghc;
    }

    if (nz > 1) {
        K = K + ghc;
    }

    return indices{I, J, K};
}
