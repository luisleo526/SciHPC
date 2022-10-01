#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

#include "scihpc/scalar_data.h"
#include "scihpc/global.h"
#include "scihpc/structured_grid.h"
#include "scihpc/vector_data.h"
#include "scihpc/runge_kutta.h"
#include "scihpc/source.h"
#include "scihpc/boundary_condition.h"
#include "scihpc/flux.h"
#include "scihpc/vtkWriter.h"

int main() {


    auto phi = scalar_data(64, 64);
    auto vel = vector_data(phi.nx, phi.ny);
    auto geo = structured_grid(axis{0.0, 1.0, phi.nx},
                               axis{0.0, 1.0, phi.ny});
    for (int i = 0; i < phi.nx; ++i) {
        for (int j = 0; j < phi.ny; ++j) {
            auto index = phi.index_mapping(i + 1, j + 1, 1);
            phi.data[index.i][index.j][index.k] = -sqrt(pow(geo.xc[i] - 0.5, 2) +
                                                        pow(geo.yc[j] - 0.5, 2)) + 0.5;
        }
    }

    auto vtk = vtkWriter(&geo, "test");
    vtk.create(0);
    vtk.add_scalar_data(&phi, "phi");
    vtk.add_vector_data(&vel, "vel");
    vtk.close();




}