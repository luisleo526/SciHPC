//
// Created by leo on 10/7/22.
//

#ifndef SCIHPC_BC_FACTORY_H
#define SCIHPC_BC_FACTORY_H

#include "global.h"
#include "vector_data.h"
#include "no_slip.h"
#include "slip.h"
#include <string>
#include "dirichlet.h"
#include "neumann.h"
#include "structured_grid.h"
#include "periodic.h"

const std::string NO_SLIP = "no_slip";
const std::string SLIP = "slip";
const std::string PERIODIC = "periodic";
const std::string NEUMANN = "neumann";
const std::string DIRICHLET = "dirichlet";

struct bc_info {
    std::string type = "NoSlip";
    DataType value = 0.0;
};

// -------------------------- BC Factory -------------------------- //
// ___face___ or ___node___ depends on the position of value.
// Example:
// x staggered -> ___face__x + ___node__y + ___node__z
// y staggered -> ___node__x + ___face__y + ___node__z
// z staggered -> ___node__x + ___node__y + ___face__z
// ---------------------------------------------------------------- //

class bc_factory {
private:
    static void check_bc(const bc_info& bc_r, const bc_info& bc_l);

public:
    bc_factory(structured_grid *geo, const bc_info &xl_bc, const bc_info &xr_bc);

    bc_factory(structured_grid *geo, const bc_info &xl_bc, const bc_info &xr_bc, const bc_info &yl_bc,
               const bc_info &yr_bc);

    bc_factory(structured_grid *geo, const bc_info &xl_bc, const bc_info &xr_bc, const bc_info &yl_bc,
               const bc_info &yr_bc, const bc_info &zl_bc, const bc_info &zr_bc);

    bc_info xlbc, xrbc, ylbc, yrbc, zlbc, zrbc;
    DataType dx, dy, dz;

    void apply_node_centered(scalar_data* f);
    void apply_node_staggered_x(scalar_data* f);
    void apply_node_staggered_y(scalar_data* f);
    void apply_node_staggered_z(scalar_data* f);


};


#endif //SCIHPC_BC_FACTORY_H
