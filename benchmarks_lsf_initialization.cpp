#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <random>
#include <ctime>

#include "scihpc/wrapper.h"
#include "scihpc/vtkWriter.h"
#include "scihpc/structured_grid.h"
#include "scihpc/source.h"
#include "scihpc/runge_kutta.h"
#include "scihpc/flux.h"
#include "scihpc/simple_bc.h"
#include "scihpc/wrapper_func.h"

int main() {

    auto phi = wrapper(new scalar_data(256, 256));
    auto vel = wrapper(new vector_data(phi.scalar->nx, phi.scalar->ny));
    auto geo = structured_grid(axis{-1.1, 1.1, phi.scalar->nx},
                               axis{-1.1, 1.1, phi.scalar->ny});
    auto solver = runge_kutta(phi.scalar->Nx, phi.scalar->Ny, phi.scalar->Nz);
    auto vtk = vtkWriter(&geo, "lsf_iris");
    auto param = new problem_parameters{};
    auto deri_solvers = derivatives_solver_alloc(phi.scalar, &geo);
    auto dummy = dummy_data_alloc(phi.scalar);
    phi.link_params(param);
    phi.link_solvers(deri_solvers);
    phi.link_dummy(dummy);

    param->dt = 0.01 * geo.h;
    param->rdt = 0.25 * geo.h;
    param->ls_width = 3.0 * geo.h;

    std::default_random_engine generator(time(NULL));
    // noise
    std::uniform_real_distribution<DataType> unif(-0.00, 0.00);
    std::uniform_real_distribution<DataType> unif2(0, 0.0);

    // Initialize phi
    for (int i = 0; i < phi.scalar->nx; ++i) {
        for (int j = 0; j < phi.scalar->ny; ++j) {
            auto index = phi.scalar->index_mapping(i + 1, j + 1, 1);
            if (geo.xc[i] > -1.0 && geo.xc[i] < 1.0) {
                if (geo.yc[j] > -0.85 * cos(0.5 * pi * geo.xc[i]) - unif(generator) and
                    geo.yc[j] < 0.85 * cos(0.5 * pi * geo.xc[i]) + unif(generator)) {
                    if (pow(geo.xc[i], 2) + pow(geo.yc[j], 2) < 0.25 * (1.0 + unif(generator))) {
                        phi.scalar->data[index.i][index.j][index.k] = 255.0 - unif2(generator);
                    } else {
                        phi.scalar->data[index.i][index.j][index.k] = -255.0 + unif2(generator);
                    }
                } else {
                    phi.scalar->data[index.i][index.j][index.k] = 255.0 - unif2(generator);
                }
            } else {
                phi.scalar->data[index.i][index.j][index.k] = 255.0 + unif2(generator);
            }

            phi.scalar->data[index.i][index.j][index.k] /= 255.0;
        }
    }
    zero_order_extrapolation(phi.scalar);

    vtk.create(0);
    vtk.add_scalar(phi.scalar, "phi");
    vtk.close();

    find_sign(&phi);
    stabilized_upon_gradient(&phi, &geo);
    int step = 0;
    do {
        store_tmp(&phi);
        if (step % 100 == 0 and step < 500) {
            find_sign(&phi);
            stabilized_upon_gradient(&phi, &geo);
        }
        solver.tvd_rk3(&phi, &vel, &geo, &identity_flux, &zero_order_extrapolation, &lsf_redistance_lambda);
        std::cout << l2norm(&phi) << std::endl;
    } while (++step < 1500 and l2norm(&phi) > 1e-6);

    vtk.create(1);
    vtk.add_scalar(phi.scalar, "phi");
    vtk.close();

    return 0;
}
