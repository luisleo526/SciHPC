//
// Created by leo on 10/7/22.
//

#include "velocities_bc.h"
#include <cassert>

void velocities_bc::check_bc(bc_info bc_r, bc_info bc_l) {

    assert(bc_r.type == PERIODIC || bc_r.type == DIRICHLET || bc_r.type == NEUMANN ||
           bc_r.type == SLIP || bc_r.type == NO_SLIP);

    assert(bc_l.type == PERIODIC || bc_l.type == DIRICHLET || bc_l.type == NEUMANN ||
           bc_l.type == SLIP || bc_l.type == NO_SLIP);

    // Check periodic
    if (bc_l.type == PERIODIC || bc_r.type == PERIODIC) {
        assert(bc_l.type == PERIODIC && bc_r.type == PERIODIC);
    }

    // Check type need value
    if (bc_l.type == DIRICHLET || bc_l.type == NEUMANN) {
        assert(bc_l.value != 0.0);
    }

    if (bc_r.type == DIRICHLET || bc_r.type == NEUMANN) {
        assert(bc_r.value != 0.0);
    }

}

velocities_bc::velocities_bc(const bc_info &xl_bc, const bc_info &xr_bc) {
    check_bc(xl_bc, xr_bc);
    xlbc = xl_bc;
    xrbc = xr_bc;
}

velocities_bc::velocities_bc(const bc_info &xl_bc, const bc_info &xr_bc, const bc_info &yl_bc, const bc_info &yr_bc) {
    check_bc(xl_bc, xr_bc);
    check_bc(yl_bc, yr_bc);
    xlbc = xl_bc;
    xrbc = xr_bc;
    ylbc = yl_bc;
    yrbc = yr_bc;
}

velocities_bc::velocities_bc(const bc_info &xl_bc, const bc_info &xr_bc, const bc_info &yl_bc, const bc_info &yr_bc,
                             const bc_info &zl_bc, const bc_info &zr_bc) {
    check_bc(xl_bc, xr_bc);
    check_bc(yl_bc, yr_bc);
    check_bc(zl_bc, zr_bc);
    xlbc = xl_bc;
    xrbc = xr_bc;
    ylbc = yl_bc;
    yrbc = yr_bc;
    zlbc = zl_bc;
    zrbc = zr_bc;
}

void velocities_bc::apply_vel_bc(vector_data *vel) const {

    if (xlbc.type == NO_SLIP) {
        no_slip_face_xl(vel);
    } else if (xlbc.type == SLIP) {
        slip_face_xl(vel);
    }

    if (xrbc.type == NO_SLIP) {
        no_slip_face_xr(vel);
    } else if (xrbc.type == SLIP) {
        slip_face_xr(vel);
    }

    if (ylbc.type == NO_SLIP) {
        no_slip_face_yl(vel);
    } else if (ylbc.type == SLIP) {
        slip_face_yl(vel);
    }

    if (yrbc.type == NO_SLIP) {
        no_slip_face_yr(vel);
    } else if (yrbc.type == SLIP) {
        slip_face_yr(vel);
    }

    if (zlbc.type == NO_SLIP) {
        no_slip_face_zl(vel);
    } else if (zlbc.type == SLIP) {
        slip_face_zl(vel);
    }

    if (zrbc.type == NO_SLIP) {
        no_slip_face_zr(vel);
    } else if (zrbc.type == SLIP) {
        slip_face_zr(vel);
    }
}

void velocities_bc::apply_nvel_bc(vector_data *nvel) const {

    if (xlbc.type == NO_SLIP) {
        no_slip_node_xl(nvel);
    } else if (xlbc.type == SLIP) {
        slip_node_xl(nvel);
    }

    if (xrbc.type == NO_SLIP) {
        no_slip_node_xr(nvel);
    } else if (xrbc.type == SLIP) {
        slip_node_xr(nvel);
    }

    if (ylbc.type == NO_SLIP) {
        no_slip_node_yl(nvel);
    } else if (ylbc.type == SLIP) {
        slip_node_yl(nvel);
    }

    if (yrbc.type == NO_SLIP) {
        no_slip_node_yr(nvel);
    } else if (yrbc.type == SLIP) {
        slip_node_yr(nvel);
    }

    if (zlbc.type == NO_SLIP) {
        no_slip_node_zl(nvel);
    } else if (zlbc.type == SLIP) {
        slip_node_zl(nvel);
    }

    if (zrbc.type == NO_SLIP) {
        no_slip_node_zr(nvel);
    } else if (zrbc.type == SLIP) {
        slip_node_zr(nvel);
    }

}

void velocities_bc::apply_velx_bc(vector_data *vel) const {


    if (xlbc.type == NO_SLIP) {
        no_slip_face_xl(&vel->x);
        if (vel->x.ndim > 1) {
            no_slip_face_xl(&vel->y);
        }
        if (vel->x.ndim > 2) {
            no_slip_face_xl(&vel->z);
        }
    } else if (xlbc.type == SLIP) {
        slip_face_xl(&vel->x);
        if (vel->x.ndim > 1) {
            slip_face_xl(&vel->y);
        }
        if (vel->x.ndim > 2) {
            slip_face_xl(&vel->z);
        }
    }

    if (xrbc.type == NO_SLIP) {
        no_slip_face_xr(&vel->x);
        if (vel->x.ndim > 1) {
            no_slip_face_xr(&vel->y);
        }
        if (vel->x.ndim > 2) {
            no_slip_face_xr(&vel->z);
        }
    } else if (xrbc.type == SLIP) {
        slip_face_xr(&vel->x);
        if (vel->x.ndim > 1) {
            slip_face_xr(&vel->y);
        }
        if (vel->x.ndim > 2) {
            slip_face_xr(&vel->z);
        }
    }

    if (ylbc.type == NO_SLIP) {
        no_slip_node_yl(vel);
    } else if (ylbc.type == SLIP) {
        slip_node_yl(vel);
    }

    if (yrbc.type == NO_SLIP) {
        no_slip_node_yr(vel);
    } else if (yrbc.type == SLIP) {
        slip_node_yr(vel);
    }

    if (zlbc.type == NO_SLIP) {
        no_slip_node_zl(vel);
    } else if (zlbc.type == SLIP) {
        slip_node_zl(vel);
    }

    if (zrbc.type == NO_SLIP) {
        no_slip_node_zr(vel);
    } else if (zrbc.type == SLIP) {
        slip_node_zr(vel);
    }

}

void velocities_bc::apply_vely_bc(vector_data *vel) const {

    if (xlbc.type == NO_SLIP) {
        no_slip_node_xl(vel);
    } else if (xlbc.type == SLIP) {
        slip_node_xl(vel);
    }

    if (xrbc.type == NO_SLIP) {
        no_slip_node_xr(vel);
    } else if (xrbc.type == SLIP) {
        slip_node_xr(vel);
    }

    if (ylbc.type == NO_SLIP) {
        if (vel->x.ndim > 1) {
            no_slip_face_yl(&vel->x);
            no_slip_face_yl(&vel->y);
        }
        if (vel->x.ndim > 2) {
            no_slip_face_yl(&vel->z);
        }
    } else if (ylbc.type == SLIP) {
        if (vel->x.ndim > 1) {
            slip_face_yl(&vel->x);
            slip_face_yl(&vel->y);
        }
        if (vel->x.ndim > 2) {
            slip_face_yl(&vel->z);
        }
    }

    if (yrbc.type == NO_SLIP) {
        if (vel->x.ndim > 1) {
            no_slip_face_yr(&vel->x);
            no_slip_face_yr(&vel->y);
        }
        if (vel->x.ndim > 2) {
            no_slip_face_yr(&vel->z);
        }
    } else if (yrbc.type == SLIP) {
        if (vel->x.ndim > 1) {
            slip_face_yr(&vel->x);
            slip_face_yr(&vel->y);
        }
        if (vel->x.ndim > 2) {
            slip_face_yr(&vel->z);
        }
    }

    if (zlbc.type == NO_SLIP) {
        no_slip_node_zl(vel);
    } else if (zlbc.type == SLIP) {
        slip_node_zl(vel);
    }

    if (zrbc.type == NO_SLIP) {
        no_slip_node_zr(vel);
    } else if (zrbc.type == SLIP) {
        slip_node_zr(vel);
    }

}

void velocities_bc::apply_velz_bc(vector_data *vel) const {
    if (xlbc.type == NO_SLIP) {
        no_slip_node_xl(vel);
    } else if (xlbc.type == SLIP) {
        slip_node_xl(vel);
    }

    if (xrbc.type == NO_SLIP) {
        no_slip_node_xr(vel);
    } else if (xrbc.type == SLIP) {
        slip_node_xr(vel);
    }

    if (ylbc.type == NO_SLIP) {
        no_slip_node_yl(vel);
    } else if (ylbc.type == SLIP) {
        slip_node_yl(vel);
    }

    if (yrbc.type == NO_SLIP) {
        no_slip_node_yr(vel);
    } else if (yrbc.type == SLIP) {
        slip_node_yr(vel);
    }

    if (zlbc.type == NO_SLIP) {
        if (vel->x.ndim > 2) {
            no_slip_face_zl(&vel->x);
            no_slip_face_zl(&vel->y);
            no_slip_face_zl(&vel->z);
        }
    } else if (zlbc.type == SLIP) {
        if (vel->x.ndim > 2) {
            slip_face_zl(&vel->x);
            slip_face_zl(&vel->y);
            slip_face_zl(&vel->z);
        }
    }

    if (zrbc.type == NO_SLIP) {
        if (vel->x.ndim > 2) {
            no_slip_face_zr(&vel->x);
            no_slip_face_zr(&vel->y);
            no_slip_face_zr(&vel->z);
        }
    } else if (zrbc.type == SLIP) {
        if (vel->x.ndim > 2) {
            slip_face_zr(&vel->x);
            slip_face_zr(&vel->y);
            slip_face_zr(&vel->z);
        }
    }

}




