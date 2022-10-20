//
// Created by 溫晧良 on 2022/10/14.
//

#include "multigrid.h"
#include <iostream>

int multigrid::check2d(std::vector<int> &nx, std::vector<int> &ny) {
    if (nx.back() % 2 == 0 && ny.back() % 2 == 0 && nx.back() / 2 >= 2 && ny.back() / 2 >= 2) {
        nx.push_back(nx.back() / 2);
        ny.push_back(ny.back() / 2);
        return 2;
    } else if (nx.back() % 3 == 0 && ny.back() % 3 == 0 && nx.back() / 3 >= 2 && ny.back() / 3 >= 2) {
        nx.push_back(nx.back() / 3);
        ny.push_back(ny.back() / 3);
        return 3;
    } else {
        return 1;
    }
}

int multigrid::check3d(std::vector<int> &nx, std::vector<int> &ny, std::vector<int> &nz) {
    if (nx.back() % 2 == 0 && ny.back() % 2 == 0 && nz.back() % 2 == 0 && nx.back() / 2 >= 2 && ny.back() / 2 >= 2 &&
        nz.back() / 2 >= 2) {
        nx.push_back(nx.back() / 2);
        ny.push_back(ny.back() / 2);
        nz.push_back(nz.back() / 2);
        return 2;
    } else if (nx.back() % 3 == 0 && ny.back() % 3 == 0 && nz.back() % 3 == 0 && nx.back() / 3 >= 2 &&
               ny.back() / 3 >= 2 && nz.back() / 3 >= 2) {
        nx.push_back(nx.back() / 3);
        ny.push_back(ny.back() / 3);
        nz.push_back(nz.back() / 3);
        return 3;
    } else {
        return 1;
    }
}

multigrid::multigrid(wrapper *var) {

    std::vector<int> degree;
    std::vector<int> nx;
    std::vector<int> ny;
    std::vector<int> nz;

    nx.push_back(var->scalar->nx);
    ny.push_back(var->scalar->ny);
    nz.push_back(var->scalar->nz);
    degree.push_back(1);

    int next_degree;
    while (true) {
        if (var->scalar->ndim == 2) {
            next_degree = check2d(nx, ny);
            if (next_degree == 1) {
                break;
            } else {
                degree.push_back(next_degree * degree.back());
            }
        } else {
            next_degree = check3d(nx, ny, nz);
            if (next_degree == 1) {
                break;
            } else {
                degree.push_back(next_degree * degree.back());
            }
        }
    }

    level_num = degree.size();
    at = new multigrid_base *[level_num];

    for (int i = 0; i < level_num; i++) {
        if (var->scalar->ndim == 2) {
            at[i] = new multigrid_base(nx[i], ny[i], degree[i], var->geo->dx * degree[i], var->geo->dy * degree[i]);
        } else {
            at[i] = new multigrid_base(nx[i], ny[i], nz[i], degree[i], var->geo->dx * degree[i],
                                       var->geo->dy * degree[i],
                                       var->geo->dz * degree[i]);
        }
    }
}

void multigrid::v_cycle() {

    for (int i = 0; i < level_num - 1; ++i) {
        at[i]->restriction(*at[i + 1]);
    }

    at[level_num - 1]->solve(1e-16);

    for (int i = int(level_num) - 1; i > 0; --i) {
        at[i]->prolongation(*at[i - 1]);
    }

}

void multigrid::full_cycle() {

    for (int i = 0; i < level_num - 1; ++i) {
        at[i]->restriction(*at[i + 1]);
    }

    at[level_num - 1]->solve(1e-16);

    for (int iter = 1; iter < level_num; ++iter) {
        for (int i = int(level_num) - 1; i < int(level_num) - 1 - iter; --i) {
            at[i]->prolongation(*at[i - 1]);
        }
        for (int i = int(level_num) - 1 - iter; i < level_num - 1; ++i) {
            at[i]->restriction(*at[i + 1]);
        }
        at[level_num - 1]->solve(1e-16);
    }

    for (int i = int(level_num) - 1; i > 0; --i) {
        at[i]->prolongation(*at[i - 1]);
    }

}
