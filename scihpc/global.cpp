//
// Created by 溫晧良 on 2022/9/28.
//

#include "global.h"

DataType ***init_array(int Nx, int Ny, int Nz) {
    auto ***arr = new DataType **[Nx];
    for (int i = 0; i < Nx; ++i) {
        arr[i] = new DataType *[Ny];
        for (int j = 0; j < Ny; ++j) {
            arr[i][j] = new DataType[Nz];
            for (int k = 0; k < Nz; ++k) {
                arr[i][j][k] = 0.0;
            }
        }
    }
    return arr;
}

void delete3d(DataType ***arr, int Nx, int Ny) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }
    delete[] arr;
}

problem_parameters *
set_parameters(DataType rho1, DataType rho2, DataType mu1, DataType mu2, DataType surface_tension, DataType gravity,
               DataType length) {
    auto param = new problem_parameters();
    param->density_ratio = rho2 / rho1;
    param->viscosity_ratio = mu2 / mu1;

    auto Lc = length;
    auto Uc = sqrt(gravity * Lc);

    param->Lc = length;
    param->Uc = sqrt(gravity * param->Lc);
    param->Tc = param->Lc / param->Uc;

    param->Reynolds_number = rho1 * param->Uc * param->Lc / mu1;
    param->Froude_number = param->Uc / sqrt(gravity * param->Lc);
    param->Weber_number = rho1 * param->Uc * param->Uc * param->Lc / surface_tension;

    return param;
}

problem_parameters *
set_parameters(DataType rho1, DataType rho2, DataType mu1, DataType mu2, DataType surface_tension, DataType gravity,
               DataType length, DataType velocity) {
    auto param = new problem_parameters();
    param->density_ratio = rho2 / rho1;
    param->viscosity_ratio = mu2 / mu1;

    auto Lc = length;
    auto Uc = sqrt(gravity * Lc);

    param->Lc = length;
    param->Uc = velocity;
    param->Tc = param->Lc / param->Uc;

    param->Reynolds_number = rho1 * param->Uc * param->Lc / mu1;
    param->Froude_number = param->Uc / sqrt(gravity * param->Lc);
    param->Weber_number = rho1 * param->Uc * param->Uc * param->Lc / surface_tension;

    return param;
}

problem_parameters *set_air_water(DataType length) {
    return set_parameters(1000.0, 1.226, 1.137e-3, 1.78e-5,
                          0.0728, 9.81, length);
}

problem_parameters *set_air_water(DataType length, DataType velocity) {
    return set_parameters(1000.0, 1.226, 1.137e-3, 1.78e-5,
                          0.0728, 9.81, length, velocity);
}
