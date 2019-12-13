//
// Created by Bowen Zhao on 9/26/19.
//
#pragma once

#include <cmath>

namespace AC{
    using Type_ValReal = double;

    constexpr Type_ValReal Val_PI = 3.14159265358979323846;

    struct type_ParaHamil{
        Type_ValReal Val_Omega_Lower;
        Type_ValReal Val_Omega_Upper;
        Type_ValReal Num_OmegaGrid;
        Type_ValReal Num_DeltaFunc;
    };


}
