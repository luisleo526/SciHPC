//
// Created by leo on 10/7/22.
//

#ifndef SCIHPC_VELOCITIES_BC_H
#define SCIHPC_VELOCITIES_BC_H

#include "global.h"
#include "vector_data.h"
#include "no_slip.h"
#include "slip.h"
#include <string>

const std::string NO_SLIP = "no_slip";
const std::string SLIP = "slip";
const std::string PERIODIC = "periodic";
const std::string NEUMANN = "neumann";
const std::string DIRICHLET = "dirichlet";

struct bc_info{
    std::string type = "NoSlip";
    DataType value = 0.0;
};

class velocities_bc {
private:
    static void check_bc(bc_info bc_r, bc_info bc_l);
public:
    velocities_bc(const bc_info& xl_bc, const bc_info& xr_bc);
    velocities_bc(const bc_info& xl_bc, const bc_info& xr_bc, const bc_info& yl_bc, const bc_info& yr_bc);
    velocities_bc(const bc_info& xl_bc, const bc_info& xr_bc, const bc_info& yl_bc, const bc_info& yr_bc, const bc_info& zl_bc, const bc_info& zr_bc);
    bc_info xlbc, xrbc, ylbc, yrbc, zlbc, zrbc;
    void apply_vel_bc(vector_data* vel) const;
    void apply_nvel_bc(vector_data* nvel) const;
    void apply_velx_bc(vector_data* vel) const;
    void apply_vely_bc(vector_data* vel) const;
    void apply_velz_bc(vector_data* vel) const;
};


#endif //SCIHPC_VELOCITIES_BC_H
