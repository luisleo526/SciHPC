//
// Created by 溫晧良 on 2022/10/1.
//

#ifndef SCIHPC_FLUX_H
#define SCIHPC_FLUX_H

#include "global.h"
#include "scalar_data.h"
#include "vector_data.h"

void burgers_flux(scalar_data *f, vector_data *vel);

void identity_flux(scalar_data *f, vector_data *vel);
void identity_flux(scalar_data *f);
void identity_flux(vector_data *vel);

void identity_old_flux(scalar_data *f, vector_data *vel);
void identity_old_flux(scalar_data *f);
void identity_old_flux(vector_data *vel);

void guess_from_history(scalar_data *f);

void identity_with_extrapolation(scalar_data *f, vector_data *vel);
void identity_with_extrapolation(scalar_data *f);
void identity_with_extrapolation(vector_data *vel);
void identity_with_extrapolation_face_x(scalar_data *f);
void identity_with_extrapolation_face_x(vector_data *f);
void identity_with_extrapolation_face_y(scalar_data *f);
void identity_with_extrapolation_face_y(vector_data *f);
void identity_with_extrapolation_face_z(scalar_data *f);
void identity_with_extrapolation_face_z(vector_data *f);
void identity_with_extrapolation_face_x(scalar_data *f, vector_data *vel);
void identity_with_extrapolation_face_y(scalar_data *f, vector_data *vel);
void identity_with_extrapolation_face_z(scalar_data *f, vector_data *vel);


void identity_old_with_extrapolation(scalar_data *f, vector_data *vel);
void identity_old_with_extrapolation(scalar_data *f);
void identity_old_with_extrapolation(vector_data *vel);
void identity_old_with_extrapolation_face_x(scalar_data *f);
void identity_old_with_extrapolation_face_x(vector_data *f);
void identity_old_with_extrapolation_face_y(scalar_data *f);
void identity_old_with_extrapolation_face_y(vector_data *f);
void identity_old_with_extrapolation_face_z(scalar_data *f);
void identity_old_with_extrapolation_face_z(vector_data *f);
void identity_old_with_extrapolation_face_x(scalar_data *f, vector_data *vel);
void identity_old_with_extrapolation_face_y(scalar_data *f, vector_data *vel);
void identity_old_with_extrapolation_face_z(scalar_data *f, vector_data *vel);

#endif //SCIHPC_FLUX_H
