//
// Created by 溫晧良 on 2022/10/14.
//

#ifndef SCIHPC_MULTIGRID_H
#define SCIHPC_MULTIGRID_H

#include <vector>
#include "wrapper.h"
#include "multigrid_base.h"

class multigrid {
private:
    int check2d(std::vector<int> &nx, std::vector<int> &ny);

    int check3d(std::vector<int> &nx, std::vector<int> &ny, std::vector<int> &nz);

public:
    multigrid(wrapper *var);
    unsigned int level_num;
    std::vector<multigrid_base> at;
    void v_cycle();
    void full_cycle();
};


#endif //SCIHPC_MULTIGRID_H
