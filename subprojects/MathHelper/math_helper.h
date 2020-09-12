#pragma once

#include <graphene.h>

void math_helper_perspective_3d (graphene_point3d_t * p1,
                              graphene_point3d_t * p2,
                              graphene_point3d_t * p3,
                              graphene_point3d_t * p4,
                              graphene_point3d_t * q1,
                              graphene_point3d_t * q2,
                              graphene_point3d_t * q3,
                              graphene_point3d_t * q4,
                              graphene_matrix_t * m);

int math_helper_singular_value_decomposition (double * A,
                                              int nrows,
                                              int ncols,
                                              double * U,
                                              double * S,
                                              double * V);

void math_helper_singular_value_decomposition_solve (double * U,
                                                     double * S,
                                                     double * V,
                                                     int nrows,
                                                     int ncols,
                                                     double * B,
                                                     double * x);

