//
// Created by 溫晧良 on 2022/10/10.
//

#ifndef SCIHPC_PERIODIC_H
#define SCIHPC_PERIODIC_H

#include "scalar_data.h"

void periodic_node_xr(scalar_data *f);
void periodic_node_xl(scalar_data *f);
void periodic_node_yr(scalar_data *f);
void periodic_node_yl(scalar_data *f);
void periodic_node_zr(scalar_data *f);
void periodic_node_zl(scalar_data *f);

void periodic_face_xr(scalar_data *f);
void periodic_face_xl(scalar_data *f);
void periodic_face_yr(scalar_data *f);
void periodic_face_yl(scalar_data *f);
void periodic_face_zr(scalar_data *f);
void periodic_face_zl(scalar_data *f);

#endif //SCIHPC_PERIODIC_H
