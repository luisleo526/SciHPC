//
// Created by leo on 10/7/22.
//

#ifndef SCIHPC_SLIP_H
#define SCIHPC_SLIP_H

#include "vector_data.h"

void slip_node_xr(scalar_data *vel);
void slip_node_xl(scalar_data *vel);
void slip_node_yr(scalar_data *vel);
void slip_node_yl(scalar_data *vel);
void slip_node_zr(scalar_data *vel);
void slip_node_zl(scalar_data *vel);

void slip_face_xr(scalar_data *u);
void slip_face_xl(scalar_data *u);
void slip_face_yr(scalar_data *v);
void slip_face_yl(scalar_data *v);
void slip_face_zr(scalar_data *w);
void slip_face_zl(scalar_data *w);


#endif //SCIHPC_SLIP_H
