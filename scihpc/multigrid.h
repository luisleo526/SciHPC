//
// Created by 溫晧良 on 2022/10/14.
//

#ifndef SCIHPC_MULTIGRID_H
#define SCIHPC_MULTIGRID_H

#include <vector>
#include "scalar_data.h"
#include "structured_grid.h"
#include "multigrid_base.h"

class multigrid {
private:
    int check2d(std::vector<int> &nx, std::vector<int> &ny);
    int check3d(std::vector<int> &nx, std::vector<int> &ny, std::vector<int> &nz);

public:
    multigrid(scalar_data *var, structured_grid* geo);

    unsigned int level_num;
    multigrid_base **at;
    void v_cycle();
    void full_cycle();
};


#endif //SCIHPC_MULTIGRID_H
