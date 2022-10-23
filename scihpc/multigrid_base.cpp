//
// Created by 溫晧良 on 2022/10/14.
//

#include "multigrid_base.h"
#include <iostream>
#include <chrono>
#include <thread>

int multigrid_base::of(int i, int j, int k) {
    return i + j * nx + k * nx * ny;
}

int multigrid_base::of(int i, int j) {
    return i + j * nx;
}

void multigrid_base::relax(int iter) {

    DataType rdot, alpha;

    // Conjugate gradient
    for (int cnt = 0; cnt < iter; ++cnt) {
        if (residual() < 1e-16) {
            break;
        }
        rdot = 0.0;
#pragma omp parallel for default(none) reduction(+:rdot)
        for (int i = 0; i < n; ++i) {
            rdot += res[i] * res[i];
        }
        Ax(res);
        alpha = 0.0;
#pragma omp parallel for default(none) reduction(+:alpha)
        for (int i = 0; i < n; ++i) {
            alpha += res[i] * buffer[i];
        }
        alpha = rdot / alpha;
#pragma omp parallel for default(none) shared(alpha)
        for (int i = 0; i < n; ++i) {
            sol[i] += alpha * res[i];
        }
        compatibility_condition(sol);
    }


}

void multigrid_base::compatibility_condition(VectorX &f) {

    if (no_compatibility) {
        return void();
    }

    DataType sum = 0.0;

#pragma omp parallel for default(none) reduction(+:sum) shared(f)
    for (int i = 0; i < n; ++i) {
        sum += f[i];
    }

    sum = sum / n;

#pragma omp parallel for default(none) shared(sum, f)
    for (int i = 0; i < n; ++i) {
        f[i] -= sum;
    }

}

DataType multigrid_base::residual() {
    res = rhs - A * sol;
    return res.norm();
}

void multigrid_base::restriction(multigrid_base *coarse) {
    relax(base_step);
    residual();
    auto r = coarse->degree / degree;

    if (ndim == 2) {
#pragma omp parallel for default(none) shared(r, coarse) collapse(2)
        for (int i = 0; i < coarse->nx; ++i) {
            for (int j = 0; j < coarse->ny; ++j) {
                coarse->rhs[coarse->of(i, j)] = 0.0;
                for (int ii = 0; ii < r; ++ii) {
                    for (int jj = 0; jj < r; ++jj) {
                        coarse->rhs[coarse->of(i, j)] += res[of(i * r + ii, j * r + jj)];
                    }
                }
                coarse->rhs[coarse->of(i, j)] /= r * r;
                coarse->sol[coarse->of(i, j)] = 0.0;
            }
        }
    } else if (ndim == 3) {
#pragma omp parallel for default(none) shared(r, coarse) collapse(3)
        for (int i = 0; i < coarse->nx; ++i) {
            for (int j = 0; j < coarse->ny; ++j) {
                for (int k = 0; k < coarse->nz; ++k) {
                    coarse->rhs[coarse->of(i, j, k)] = 0.0;
                    for (int ii = 0; ii < r; ++ii) {
                        for (int jj = 0; jj < r; ++jj) {
                            for (int kk = 0; kk < r; ++kk) {
                                coarse->rhs[coarse->of(i, j, k)] += res[of(i * r + ii, j * r + jj, k * r + kk)];
                            }
                        }
                    }
                    coarse->rhs[coarse->of(i, j, k)] /= r * r * r;
                    coarse->sol[coarse->of(i, j, k)] = 0.0;
                }
            }
        }
    }
    coarse->compatibility_condition(coarse->rhs);
}

void multigrid_base::prolongation(multigrid_base *dense) {
    auto r = degree / dense->degree;

    DataType tmp1, tmp2;
    if (ndim == 2) {
#pragma omp parallel for default(none) shared(r, dense) collapse(2)
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {

                DataType fx, fy;

                if (i > 0 && i < nx - 1) {
                    fx = 0.5 * (sol[of(i + 1, j)] - sol[of(i - 1, j)]) / dx;
                } else {
                    if (i == 0) {
                        fx = (sol[of(i + 1, j)] - sol[of(i, j)]) / dx;
                    } else {
                        fx = (sol[of(i, j)] - sol[of(i - 1, j)]) / dx;
                    }
                }
                if (j > 0 && j < ny - 1) {
                    fy = 0.5 * (sol[of(i, j + 1)] - sol[of(i, j - 1)]) / dy;
                } else {
                    if (j == 0) {
                        fy = (sol[of(i, j + 1)] - sol[of(i, j)]) / dy;
                    } else {
                        fy = (sol[of(i, j)] - sol[of(i, j - 1)]) / dy;
                    }
                }

                for (int ii = 0; ii < r; ++ii) {
                    for (int jj = 0; jj < r; ++jj) {
                        auto x = (0.5 + ii) * dense->dx;
                        auto y = (0.5 + jj) * dense->dy;
                        tmp1 = x;
                        tmp2 = y;
                        dense->sol[dense->of(i * r + ii, j * r + jj)] +=
                                sol[of(i, j)] + fx * (x - 0.5 * dx) + fy * (y - 0.5 * dy);
                    }
                }

            }
        }
    } else if (ndim == 3) {
#pragma omp parallel for default(none) shared(r, dense) collapse(3)
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {

                    DataType fx, fy, fz;

                    if (i > 0 && i < nx - 1) {
                        fx = 0.5 * (sol[of(i + 1, j, k)] - sol[of(i - 1, j, k)]) / dx;
                    } else {
                        if (i == 0) {
                            fx = (sol[of(i + 1, j, k)] - sol[of(i, j, k)]) / dx;
                        } else {
                            fx = (sol[of(i, j, k)] - sol[of(i - 1, j, k)]) / dx;
                        }
                    }
                    if (j > 0 && j < ny - 1) {
                        fy = 0.5 * (sol[of(i, j + 1, k)] - sol[of(i, j - 1, k)]) / dy;
                    } else {
                        if (j == 0) {
                            fy = (sol[of(i, j + 1, k)] - sol[of(i, j, k)]) / dy;
                        } else {
                            fy = (sol[of(i, j, k)] - sol[of(i, j - 1, k)]) / dy;
                        }
                    }
                    if (k > 0 && k < nz - 1) {
                        fz = 0.5 * (sol[of(i, j, k + 1)] - sol[of(i, j, k - 1)]) / dz;
                    } else {
                        if (k == 0) {
                            fz = (sol[of(i, j, k + 1)] - sol[of(i, j, k)]) / dz;
                        } else {
                            fz = (sol[of(i, j, k)] - sol[of(i, j, k - 1)]) / dz;
                        }
                    }

                    for (int ii = 0; ii < r; ++ii) {
                        for (int jj = 0; jj < r; ++jj) {
                            for (int kk = 0; kk < r; ++kk) {
                                auto x = (0.5 + ii) * dense->dx;
                                auto y = (0.5 + jj) * dense->dy;
                                auto z = (0.5 + kk) * dense->dz;
                                dense->sol[dense->of(i * r + ii, j * r + jj, k * r + kk)] +=
                                        sol[of(i, j, k)] + fx * (x - 0.5 * dx) + fy * (y - 0.5 * dy) +
                                        fz * (z - 0.5 * dz);
                            }
                        }
                    }
                }
            }
        }
    }
    dense->relax(base_step);
//    std::cout << tmp1 << ", " << tmp2 << std::endl;
}

void multigrid_base::solve(DataType tol) {

    if (solver == nullptr) {
        solver = new SpSolver();
        solver->analyzePattern(A);
        solver->factorize(A);
    }

    sol = solver->solve(rhs);
    compatibility_condition(sol);
    while (residual() > 1e-10) {
        sol += solver->solve(res);
        compatibility_condition(sol);
    }

}

void multigrid_base::init_full() {

    if (ndim == 2) {
        for (int i = 0; i < n; ++i) {
            A.coeffRef(i, i) = -(2.0 / dx / dx + 2.0 / dy / dy);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            A.coeffRef(i, i) = -(2.0 / dx / dx + 2.0 / dy / dy + 2.0 / dz / dz);
        }
    }

}

multigrid_base::multigrid_base(int _nx, int _ny, int _nz, int _degree, DataType _dx, DataType _dy, DataType _dz) {

    nx = _nx;
    ny = _ny;
    nz = _nz;
    degree = _degree;
    ndim = 3;
    dx = _dx;
    dy = _dy;
    dz = _dz;
    n = nx * ny * nz;

    rhs = VectorX(n);
    sol = VectorX(n);
    res = VectorX(n);
    buffer = VectorX(n);
    A = SMatrix(n, n);
    solver = nullptr;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                DataType cc = 0.0;
                if (i > 0) {
                    A.insert(of(i, j, k), of(i - 1, j, k)) = 1.0 / dx / dx;
                    cc += 1.0 / dx / dx;
                }
                if (i < nx - 1) {
                    A.insert(of(i, j, k), of(i + 1, j, k)) = 1.0 / dx / dx;
                    cc += 1.0 / dx / dx;
                }
                if (j > 0) {
                    A.insert(of(i, j, k), of(i, j - 1, k)) = 1.0 / dy / dy;
                    cc += 1.0 / dy / dy;
                }
                if (j < ny - 1) {
                    A.insert(of(i, j, k), of(i, j + 1, k)) = 1.0 / dy / dy;
                    cc += 1.0 / dy / dy;
                }
                if (k > 0) {
                    A.insert(of(i, j, k), of(i, j, k - 1)) = 1.0 / dz / dz;
                    cc += 1.0 / dz / dz;
                }
                if (k < nz - 1) {
                    A.insert(of(i, j, k), of(i, j, k + 1)) = 1.0 / dz / dz;
                    cc += 1.0 / dz / dz;
                }
                A.insert(of(i, j, k), of(i, j, k)) = -cc;
            }
        }
    }

    A.makeCompressed();
}

multigrid_base::multigrid_base(int _nx, int _ny, int _degree, DataType _dx, DataType _dy) {
    nx = _nx;
    ny = _ny;
    nz = 1;
    degree = _degree;
    ndim = 2;
    dx = _dx;
    dy = _dy;
    dz = 1.0;
    n = nx * ny;

    rhs = VectorX(n);
    sol = VectorX(n);
    res = VectorX(n);
    buffer = VectorX(n);
    A = SMatrix(n, n);
    solver = nullptr;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            DataType cc = 0.0;
            if (i > 0) {
                A.insert(of(i, j), of(i - 1, j)) = 1.0 / dx / dx;
                cc += 1.0 / dx / dx;
            }
            if (i < nx - 1) {
                A.insert(of(i, j), of(i + 1, j)) = 1.0 / dx / dx;
                cc += 1.0 / dx / dx;
            }
            if (j > 0) {
                A.insert(of(i, j), of(i, j - 1)) = 1.0 / dy / dy;
                cc += 1.0 / dy / dy;
            }
            if (j < ny - 1) {
                A.insert(of(i, j), of(i, j + 1)) = 1.0 / dy / dy;
                cc += 1.0 / dy / dy;
            }
            A.insert(of(i, j), of(i, j)) = -cc;
        }
    }
    A.makeCompressed();
}

void multigrid_base::Ax(const VectorX &x) {

    buffer = A * x;

}

