//
// Created by leo on 10/7/22.
//

#include "bc_factory.h"
#include <cassert>

void bc_factory::check_bc(const bc_info &bc_r, const bc_info &bc_l) {

    assert(bc_r.type == PERIODIC || bc_r.type == DIRICHLET || bc_r.type == NEUMANN ||
           bc_r.type == SLIP || bc_r.type == NO_SLIP);

    assert(bc_l.type == PERIODIC || bc_l.type == DIRICHLET || bc_l.type == NEUMANN ||
           bc_l.type == SLIP || bc_l.type == NO_SLIP);

    // Check periodic
    if (bc_l.type == PERIODIC || bc_r.type == PERIODIC) {
        assert(bc_l.type == PERIODIC && bc_r.type == PERIODIC);
    }

}

bc_factory::bc_factory(structured_grid *geo, const bc_info &xl_bc, const bc_info &xr_bc) {
    check_bc(xl_bc, xr_bc);
    xlbc = xl_bc;
    xrbc = xr_bc;
    dx = geo->dx;
    dy = geo->dy;
    dz = geo->dz;
}

bc_factory::bc_factory(structured_grid *geo, const bc_info &xl_bc, const bc_info &xr_bc, const bc_info &yl_bc,
                       const bc_info &yr_bc) {
    check_bc(xl_bc, xr_bc);
    check_bc(yl_bc, yr_bc);
    xlbc = xl_bc;
    xrbc = xr_bc;
    ylbc = yl_bc;
    yrbc = yr_bc;
    dx = geo->dx;
    dy = geo->dy;
    dz = geo->dz;
}

bc_factory::bc_factory(structured_grid *geo, const bc_info &xl_bc, const bc_info &xr_bc, const bc_info &yl_bc,
                       const bc_info &yr_bc, const bc_info &zl_bc, const bc_info &zr_bc) {
    check_bc(xl_bc, xr_bc);
    check_bc(yl_bc, yr_bc);
    check_bc(zl_bc, zr_bc);
    xlbc = xl_bc;
    xrbc = xr_bc;
    ylbc = yl_bc;
    yrbc = yr_bc;
    zlbc = zl_bc;
    zrbc = zr_bc;
    dx = geo->dx;
    dy = geo->dy;
    dz = geo->dz;
}

void bc_factory::apply_node_centered(scalar_data *f) {

    if (xrbc.type == NO_SLIP) {
        no_slip_node_xr(f);
    } else if (xrbc.type == SLIP) {
        slip_node_xr(f);
    } else if (xrbc.type == DIRICHLET) {
        Dirichlet_node_xr(f, xrbc.value);
    } else if (xrbc.type == NEUMANN) {
        Neumann_node_xr(f, xrbc.value, dx);
    } else if (xrbc.type == PERIODIC){
        periodic_node_xr(f);
    }

    if (xlbc.type == NO_SLIP) {
        no_slip_node_xl(f);
    } else if (xlbc.type == SLIP) {
        slip_node_xl(f);
    } else if (xlbc.type == DIRICHLET) {
        Dirichlet_node_xl(f, xlbc.value);
    } else if (xlbc.type == NEUMANN) {
        Neumann_node_xl(f, xlbc.value, dx);
    } else if (xlbc.type == PERIODIC){
        periodic_node_xl(f);
    }

    if (f->ndim > 1) {
        if (yrbc.type == NO_SLIP) {
            no_slip_node_yr(f);
        } else if (yrbc.type == SLIP) {
            slip_node_yr(f);
        } else if (yrbc.type == DIRICHLET) {
            Dirichlet_node_yr(f, yrbc.value);
        } else if (yrbc.type == NEUMANN) {
            Neumann_node_yr(f, yrbc.value, dy);
        } else if (yrbc.type == PERIODIC){
            periodic_node_yr(f);
        }

        if (ylbc.type == NO_SLIP) {
            no_slip_node_yl(f);
        } else if (ylbc.type == SLIP) {
            slip_node_yl(f);
        } else if (ylbc.type == DIRICHLET) {
            Dirichlet_node_yl(f, ylbc.value);
        } else if (ylbc.type == NEUMANN) {
            Neumann_node_yl(f, ylbc.value, dy);
        } else if (ylbc.type == PERIODIC){
            periodic_node_yl(f);
        }
    }

    if (f->ndim > 2) {
        if (zrbc.type == NO_SLIP) {
            no_slip_node_zr(f);
        } else if (zrbc.type == SLIP) {
            slip_node_zr(f);
        } else if (zrbc.type == DIRICHLET) {
            Dirichlet_node_zr(f, zrbc.value);
        } else if (zrbc.type == NEUMANN) {
            Neumann_node_zr(f, zrbc.value, dz);
        } else if (zrbc.type == PERIODIC){
            periodic_node_zr(f);
        }

        if (zlbc.type == NO_SLIP) {
            no_slip_node_zl(f);
        } else if (zlbc.type == SLIP) {
            slip_node_zl(f);
        } else if (zlbc.type == DIRICHLET) {
            Dirichlet_node_zl(f, zlbc.value);
        } else if (zlbc.type == NEUMANN) {
            Neumann_node_zl(f, zlbc.value, dz);
        } else if (zlbc.type == PERIODIC){
            periodic_node_zl(f);
        }
    }
}

void bc_factory::apply_node_staggered_x(scalar_data *f) {

    if (xrbc.type == NO_SLIP) {
        no_slip_face_xr(f);
    } else if (xrbc.type == SLIP) {
        slip_face_xr(f);
    } else if (xrbc.type == DIRICHLET) {
        Dirichlet_face_xr(f, xrbc.value);
    } else if (xrbc.type == NEUMANN) {
        Neumann_face_xr(f, xrbc.value, dx);
    } else if (xrbc.type == PERIODIC){
        periodic_face_xr(f);
    }

    if (xlbc.type == NO_SLIP) {
        no_slip_face_xl(f);
    } else if (xlbc.type == SLIP) {
        slip_face_xl(f);
    } else if (xlbc.type == DIRICHLET) {
        Dirichlet_face_xl(f, xlbc.value);
    } else if (xlbc.type == NEUMANN) {
        Neumann_face_xl(f, xlbc.value, dx);
    } else if (xlbc.type == PERIODIC){
        periodic_face_xl(f);
    }

    if (f->ndim > 1) {
        if (yrbc.type == NO_SLIP) {
            no_slip_node_yr(f);
        } else if (yrbc.type == SLIP) {
            slip_node_yr(f);
        } else if (yrbc.type == DIRICHLET) {
            Dirichlet_node_yr(f, yrbc.value);
        } else if (yrbc.type == NEUMANN) {
            Neumann_node_yr(f, yrbc.value, dy);
        } else if (yrbc.type == PERIODIC){
            periodic_node_yr(f);
        }

        if (ylbc.type == NO_SLIP) {
            no_slip_node_yl(f);
        } else if (ylbc.type == SLIP) {
            slip_node_yl(f);
        } else if (ylbc.type == DIRICHLET) {
            Dirichlet_node_yl(f, ylbc.value);
        } else if (ylbc.type == NEUMANN) {
            Neumann_node_yl(f, ylbc.value, dy);
        } else if (ylbc.type == PERIODIC){
            periodic_node_yl(f);
        }
    }

    if (f->ndim > 2) {
        if (zrbc.type == NO_SLIP) {
            no_slip_node_zr(f);
        } else if (zrbc.type == SLIP) {
            slip_node_zr(f);
        } else if (zrbc.type == DIRICHLET) {
            Dirichlet_node_zr(f, zrbc.value);
        } else if (zrbc.type == NEUMANN) {
            Neumann_node_zr(f, zrbc.value, dz);
        } else if (zrbc.type == PERIODIC){
            periodic_node_zr(f);
        }

        if (zlbc.type == NO_SLIP) {
            no_slip_node_zl(f);
        } else if (zlbc.type == SLIP) {
            slip_node_zl(f);
        } else if (zlbc.type == DIRICHLET) {
            Dirichlet_node_zl(f, zlbc.value);
        } else if (zlbc.type == NEUMANN) {
            Neumann_node_zl(f, zlbc.value, dz);
        } else if (zlbc.type == PERIODIC){
            periodic_node_zl(f);
        }
    }
}

void bc_factory::apply_node_staggered_y(scalar_data *f) {

    if (xrbc.type == NO_SLIP) {
        no_slip_node_xr(f);
    } else if (xrbc.type == SLIP) {
        slip_node_xr(f);
    } else if (xrbc.type == DIRICHLET) {
        Dirichlet_node_xr(f, xrbc.value);
    } else if (xrbc.type == NEUMANN) {
        Neumann_node_xr(f, xrbc.value, dx);
    } else if (xrbc.type == PERIODIC){
        periodic_node_xr(f);
    }

    if (xlbc.type == NO_SLIP) {
        no_slip_node_xl(f);
    } else if (xlbc.type == SLIP) {
        slip_node_xl(f);
    } else if (xlbc.type == DIRICHLET) {
        Dirichlet_node_xl(f, xlbc.value);
    } else if (xlbc.type == NEUMANN) {
        Neumann_node_xl(f, xlbc.value, dx);
    } else if (xlbc.type == PERIODIC){
        periodic_node_xl(f);
    }

    if (f->ndim > 1) {
        if (yrbc.type == NO_SLIP) {
            no_slip_face_yr(f);
        } else if (yrbc.type == SLIP) {
            slip_face_yr(f);
        } else if (yrbc.type == DIRICHLET) {
            Dirichlet_face_yr(f, yrbc.value);
        } else if (yrbc.type == NEUMANN) {
            Neumann_face_yr(f, yrbc.value, dy);
        } else if (yrbc.type == PERIODIC){
            periodic_face_yr(f);
        }

        if (ylbc.type == NO_SLIP) {
            no_slip_face_yl(f);
        } else if (ylbc.type == SLIP) {
            slip_face_yl(f);
        } else if (ylbc.type == DIRICHLET) {
            Dirichlet_face_yl(f, ylbc.value);
        } else if (ylbc.type == NEUMANN) {
            Neumann_face_yl(f, ylbc.value, dy);
        } else if (ylbc.type == PERIODIC){
            periodic_face_yl(f);
        }
    }

    if (f->ndim > 2) {
        if (zrbc.type == NO_SLIP) {
            no_slip_node_zr(f);
        } else if (zrbc.type == SLIP) {
            slip_node_zr(f);
        } else if (zrbc.type == DIRICHLET) {
            Dirichlet_node_zr(f, zrbc.value);
        } else if (zrbc.type == NEUMANN) {
            Neumann_node_zr(f, zrbc.value, dz);
        } else if (zrbc.type == PERIODIC){
            periodic_node_zr(f);
        }

        if (zlbc.type == NO_SLIP) {
            no_slip_node_zl(f);
        } else if (zlbc.type == SLIP) {
            slip_node_zl(f);
        } else if (zlbc.type == DIRICHLET) {
            Dirichlet_node_zl(f, zlbc.value);
        } else if (zlbc.type == NEUMANN) {
            Neumann_node_zl(f, zlbc.value, dz);
        } else if (zlbc.type == PERIODIC){
            periodic_node_zl(f);
        }
    }
}

void bc_factory::apply_node_staggered_z(scalar_data *f) {

    if (xrbc.type == NO_SLIP) {
        no_slip_node_xr(f);
    } else if (xrbc.type == SLIP) {
        slip_node_xr(f);
    } else if (xrbc.type == DIRICHLET) {
        Dirichlet_node_xr(f, xrbc.value);
    } else if (xrbc.type == NEUMANN) {
        Neumann_node_xr(f, xrbc.value, dx);
    } else if (xrbc.type == PERIODIC){
        periodic_node_xr(f);
    }

    if (xlbc.type == NO_SLIP) {
        no_slip_node_xl(f);
    } else if (xlbc.type == SLIP) {
        slip_node_xl(f);
    } else if (xlbc.type == DIRICHLET) {
        Dirichlet_node_xl(f, xlbc.value);
    } else if (xlbc.type == NEUMANN) {
        Neumann_node_xl(f, xlbc.value, dx);
    } else if (xlbc.type == PERIODIC){
        periodic_node_xl(f);
    }

    if (f->ndim > 1) {
        if (yrbc.type == NO_SLIP) {
            no_slip_node_yr(f);
        } else if (yrbc.type == SLIP) {
            slip_node_yr(f);
        } else if (yrbc.type == DIRICHLET) {
            Dirichlet_node_yr(f, yrbc.value);
        } else if (yrbc.type == NEUMANN) {
            Neumann_node_yr(f, yrbc.value, dy);
        } else if (yrbc.type == PERIODIC){
            periodic_node_yr(f);
        }

        if (ylbc.type == NO_SLIP) {
            no_slip_node_yl(f);
        } else if (ylbc.type == SLIP) {
            slip_node_yl(f);
        } else if (ylbc.type == DIRICHLET) {
            Dirichlet_node_yl(f, ylbc.value);
        } else if (ylbc.type == NEUMANN) {
            Neumann_node_yl(f, ylbc.value, dy);
        } else if (ylbc.type == PERIODIC){
            periodic_node_yl(f);
        }
    }

    if (f->ndim > 2) {
        if (zrbc.type == NO_SLIP) {
            no_slip_face_zr(f);
        } else if (zrbc.type == SLIP) {
            slip_face_zr(f);
        } else if (zrbc.type == DIRICHLET) {
            Dirichlet_face_zr(f, zrbc.value);
        } else if (zrbc.type == NEUMANN) {
            Neumann_face_zr(f, zrbc.value, dz);
        } else if (zrbc.type == PERIODIC){
            periodic_face_zr(f);
        }

        if (zlbc.type == NO_SLIP) {
            no_slip_face_zl(f);
        } else if (zlbc.type == SLIP) {
            slip_face_zl(f);
        } else if (zlbc.type == DIRICHLET) {
            Dirichlet_face_zl(f, zlbc.value);
        } else if (zlbc.type == NEUMANN) {
            Neumann_face_zl(f, zlbc.value, dz);
        } else if (zlbc.type == PERIODIC){
            periodic_face_zl(f);
        }
    }
}





