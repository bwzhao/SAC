//
// Created by Bowen Zhao on 9/26/19.
//
#pragma once

#include <cmath>
#include <armadillo>

namespace AC{

    using Type_ValReal = double;

    constexpr Type_ValReal Val_PI = 3.14159265358979323846;

    constexpr double LOWER_ACCEPT = 0.4;
    constexpr double UPPER_ACCEPT = 0.5;
    constexpr double RATIO_STEPCHANGE = 1.5;
    constexpr double RATIO_TEMPCHANGE = 1.1;
    constexpr Type_ValReal THRESHOD_CHI2 = 10e-3;

    constexpr Type_ValReal START_THETA = 100.;

    constexpr int MAX_ANNEALING = 500;
    constexpr int NUM_MEASURE_ANNEALING = 20000;

    constexpr int EQ_STEPS = 1000;
    constexpr int EQ_BINS = 50;

//    constexpr Type_ValReal G_0 = 0.7957765907451;
    constexpr Type_ValReal G_0 = 0.5;

    struct type_ParaHamil{
        Type_ValReal Val_Omega_Lower;
        Type_ValReal Val_Omega_Upper;
        Type_ValReal Num_OmegaGrid;
        Type_ValReal Num_DeltaFunc;
    };

    // index:
    // 0: val_theta
    // 1: val_SpecFunc
    class Class_SpecFunc;
    using Type_AnnealData = std::tuple<Type_ValReal, Type_ValReal>;

    // Info (index, pos_old, pos_new
    using Type_UpdateInfo_Small = std::tuple<int, int, int>;

    // Info (pos_old, pos_new)
    using Type_UpdateInfo_Large = std::tuple<int, int>;

    using Type_Spectral = arma::Col<AC::Type_ValReal>;
}
