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
#include "scihpc/wrapper_func.h"

int main() {

    auto geo = structured_grid(axis{0.0, 1.0, 128},
                               axis{0.0, 1.0, 128});

    auto phi = wrapper(true, &geo,
                       bc_info{NEUMANN}, bc_info{NEUMANN},
                       bc_info{NEUMANN}, bc_info{NEUMANN});
    auto vel = wrapper(false, &geo,
                       bc_info{NEUMANN}, bc_info{NEUMANN},
                       bc_info{NEUMANN}, bc_info{NEUMANN});
    auto solver = runge_kutta(phi.scalar->Nx, phi.scalar->Ny, phi.scalar->Nz);
    auto vtk = vtkWriter(&geo, "zalesak_disk");
    auto param = new problem_parameters{};
    auto deri_solvers = SharedSolvers_alloc(phi.scalar, &geo);
    auto dummy = dummy_data_alloc(phi.scalar);
    phi.link_params(param);
    phi.link_solvers(deri_solvers);
    phi.link_dummy(dummy);

    param->dt = 0.1 * geo.h;
    param->rdt = 0.5 * geo.h;
    param->ls_width = 1.5 * geo.h;

    auto period = 2.0 * pi;


    int refine = 30;
    // Initialize phi
    for (int i = 0; i < phi.scalar->nx; ++i) {
        for (int j = 0; j < phi.scalar->ny; ++j) {
            auto index = phi.scalar->index_mapping(i + 1, j + 1, 1);

            phi.scalar->data[index.i][index.j][index.k] = 0.0;
            vel.vector->x.data[index.i][index.j][index.k] = -(geo.yc[j] - 0.5);
            vel.vector->y.data[index.i][index.j][index.k] = +(geo.xc[i] - 0.5);

            for (int ii = 0; ii < refine; ++ii) {
                for (int jj = 0; jj < refine; ++jj) {
                    auto x = geo.x[i] + ii * geo.dx / refine;
                    auto y = geo.y[j] + jj * geo.dy / refine;
                    if (pow(x - 0.5, 2) + pow(y - 0.75, 2) <= 0.15 * 0.15) {
                        if (!(fabs(x - 0.5) <= 0.025 && y <= 0.75 + 0.15 / 2.0)) {
                            phi.scalar->data[index.i][index.j][index.k] += 1.0 / refine / refine;
                        }
                    }
                }
            }

            phi.scalar->data[index.i][index.j][index.k] = 2.0 * phi.scalar->data[index.i][index.j][index.k] - 1.0;
        }
    }
    phi.apply_scalar_bc();
    vel.apply_nvel_bc();

    vtk.create(-1);
    vtk.add_scalar(phi.scalar, "phi");
    vtk.close();

    find_sign(&phi);
    stabilized_upon_gradient(&phi);
    DataType error = 0.0;
    int step = 0;
    do {
        store_tmp(&phi);
        solver.tvd_rk3(&phi, &vel, &identity_flux, &lsf_redistance_lambda);
        error = l2norm(&phi);
        std::cout << error << std::endl;
    } while (++step * param->rdt < 1.0 and error > 1e-6);


    vtk.create(0);
    vtk.add_scalar(phi.scalar, "phi");
    vtk.close();

    param->lsf_mass0 = lsf_mass(&phi);

    step = 0;
    auto instep = 0;
    auto pltid = 1;
    do {
        solver.tvd_rk3(&phi, &vel, &identity_flux, &convection);

        instep = 0;
        find_sign(&phi);
        do {
            solver.tvd_rk3(&phi, &vel, &identity_flux, &lsf_redistance_lambda);
        } while (++instep * param->rdt < 1.5 * param->ls_width);

        if (++step * param->dt >= pltid * period / 4.0) {
            vtk.create(pltid);
            vtk.add_scalar(phi.scalar, "phi");
            vtk.close();
            ++pltid;
        }
        if (step % 10 == 0) {
            std::cout << "time = " << step * param->dt
                      << ", mass loss = " << (1.0 - lsf_mass(&phi) / param->lsf_mass0) * 100 << "%"
                      << std::endl;
        }
    } while (step * phi.params->dt < period);

    return 0;
}
