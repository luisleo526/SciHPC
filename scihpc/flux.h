//
// Created by 溫晧良 on 2022/10/1.
//

#ifndef SCIHPC_FLUX_H
#define SCIHPC_FLUX_H

#include "global.h"
#include "scalar_data.h"
#include "vector_data.h"

void identity_flux(scalar_data *f, vector_data *vel);
void identity_flux(scalar_data *f);
void identity_flux(vector_data *vel);

void burgers_flux(scalar_data *f, vector_data *vel);


#endif //SCIHPC_FLUX_H
