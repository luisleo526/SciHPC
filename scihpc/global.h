//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_GLOBAL_H
#define SCIHPC_GLOBAL_H

#include <cmath>
#include <string>

typedef double DataType;

const DataType pi = acos(static_cast<DataType>(-1.0));
const DataType epsilon = 1e-8;

struct indices {
    int i, j, k;
};

struct axis {
    DataType start, end;
    int n;
};

struct problem_parameters {
    int iter = 0;
    DataType t = 0.0;
    DataType dt = 1.0;
    DataType max_CFL = 0.01;
    DataType stable_CFL = 0.001;
    DataType rdt{};

    DataType density_ratio = 1.0;
    DataType viscosity_ratio = 1.0;

    DataType ls_width{};
    DataType lsf_mass0{};
    bool positive_ref = true;

    DataType Lc{};
    DataType Uc{};
    DataType Tc{};
    DataType Reynolds_number{};
    DataType Froude_number = -1.0;
    DataType Weber_number = -1.0;

    DataType ppe_tol = 1e-4;
    DataType ppe_tol2 = 1e-7;
    int ppe_max_iter = 1000000;
    DataType ppe_omega = 1.5;
    int ppe_initer = 10;
};

DataType ***init_array(int Nx, int Ny, int Nz);

void delete3d(DataType ***arr, int Nx, int Ny);

problem_parameters *set_parameters(DataType rho1, DataType rho2,
                                   DataType mu1, DataType mu2,
                                   DataType surface_tension,
                                   DataType gravity, DataType length);

problem_parameters *set_parameters(DataType rho1, DataType rho2,
                                   DataType mu1, DataType mu2,
                                   DataType surface_tension,
                                   DataType gravity, DataType length, DataType velocity);

problem_parameters *set_air_water(DataType length);

problem_parameters *set_air_water(DataType length, DataType velocity);

#endif //SCIHPC_GLOBAL_H
