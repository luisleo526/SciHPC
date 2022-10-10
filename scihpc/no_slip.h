//
// Created by 溫晧良 on 2022/10/5.
//

#ifndef SCIHPC_NO_SLIP_H
#define SCIHPC_NO_SLIP_H

#include "vector_data.h"

void no_slip_node_xr(scalar_data *vel);
void no_slip_node_xl(scalar_data *vel);
void no_slip_node_yr(scalar_data *vel);
void no_slip_node_yl(scalar_data *vel);
void no_slip_node_zr(scalar_data *vel);
void no_slip_node_zl(scalar_data *vel);

void no_slip_face_xr(scalar_data *u);
void no_slip_face_xl(scalar_data *u);
void no_slip_face_yr(scalar_data *v);
void no_slip_face_yl(scalar_data *v);
void no_slip_face_zr(scalar_data *w);
void no_slip_face_zl(scalar_data *w);

#endif //SCIHPC_NO_SLIP_H
