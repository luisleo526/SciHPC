#include "scihpc/wrapper.h"
#include "scihpc/structured_grid.h"
#include "scihpc/multigrid.h"
#include <iostream>
#include "scihpc/vtkWriter.h"

int main() {

    omp_set_num_threads(8);
    auto geo = structured_grid(axis{0.0, 1.0, 64},
                               axis{0.0, 1.0, 64},
                               axis{0.0, 1.0, 64});
    auto phi = wrapper(true, &geo,
                       bc_info{NEUMANN}, bc_info{NEUMANN},
                       bc_info{NEUMANN}, bc_info{NEUMANN},
                       bc_info{NEUMANN}, bc_info{NEUMANN});

    auto mg = multigrid(phi.scalar, &geo);

    auto vtk = vtkWriter(&geo, "multigrid_test");

    for (int i = 0; i < mg.level_num; ++i) {
        // For Dirichlet BC, the boundary value is 0.0
        mg.at[i]->no_compatibility = true;
        mg.at[i]->init_DirichletBC();
        std::cout << mg.at[i]->ndim << " " << mg.at[i]->nx << " " << mg.at[i]->ny << " " << mg.at[i]->nz << std::endl;
    }

    std::cout << "Assigning initial values..." << std::endl;

    for (int i = 0; i < phi.scalar->nx; ++i) {
        for (int j = 0; j < phi.scalar->ny; ++j) {
            for (int k = 0; k < phi.scalar->nz; ++k) {
                mg.at[0]->rhs[mg.at[0]->of(i, j, k)] =
                        -3.0 * M_PI * M_PI * sin(M_PI * geo.xc[i]) * sin(M_PI * geo.yc[j]) * sin(M_PI * geo.zc[k]);
                auto index = phi.scalar->index_mapping(i + 1, j + 1, k + 1);
                phi.scalar->data[index.i][index.j][index.k] =
                        sin(M_PI * geo.xc[i]) * sin(M_PI * geo.yc[j]) * sin(M_PI * geo.zc[k]);
            }
        }
    }

    int step = 0;
    while (mg.at[0]->residual() > 1e-6) {
        mg.full_cycle();
        std::cout << step++ << "," << mg.at[0]->residual() << "," << mg.at[mg.level_num - 1]->residual() << std::endl;
    }

    vtk.create(0);
    vtk.add_scalar(phi.scalar, "Exact");

    DataType error = 0;
    for (int i = 0; i < phi.scalar->nx; ++i) {
        for (int j = 0; j < phi.scalar->ny; ++j) {
            for (int k = 0; k < phi.scalar->nz; ++k) {
                auto index = phi.scalar->index_mapping(i + 1, j + 1, k + 1);
                error += pow(phi.scalar->data[index.i][index.j][index.k] - mg.at[0]->sol[mg.at[0]->of(i, j)], 2);
                phi.scalar->data[index.i][index.j][index.k] = mg.at[0]->sol[mg.at[0]->of(i, j)];
            }
        }
    }

    std::cout << "Error: " << sqrt(error / phi.scalar->nx / phi.scalar->ny / phi.scalar->nz) << std::endl;

    vtk.add_scalar(phi.scalar, "Numerical");
    vtk.close();

}
