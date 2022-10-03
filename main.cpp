#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <omp.h>

#include "scihpc/wrapper.h"
#include "scihpc/vtkWriter.h"
#include "scihpc/structured_grid.h"
#include "scihpc/source.h"
#include "scihpc/runge_kutta.h"

int main() {

    auto phi = wrapper(new scalar_data(64, 64));
    auto vel = wrapper(new vector_data(phi.scalar->nx, phi.scalar->ny));
    auto geo = structured_grid(axis{0.0, 1.0, phi.scalar->nx},
                               axis{0.0, 1.0, phi.scalar->ny});
    auto solver = runge_kutta(phi.scalar->Nx, phi.scalar->Ny, phi.scalar->Nz);
    auto vtk = vtkWriter(&geo, "lsf_initialization");
    auto param = new problem_parameters{};
    auto deri_solvers = derivatives_solver_alloc(phi.scalar, &geo);
    auto dummy = dummy_data_alloc(phi.scalar);
    phi.link_params(param);
    phi.link_solvers(deri_solvers);
    phi.link_dummy(dummy);

    // Initialize phi
    for (int i = 0; i < phi.scalar->nx; ++i) {
        for (int j = 0; j < phi.scalar->ny; ++j) {
            
        }
    }


    std::cout << phi.params->dt;

    return 0;
}
