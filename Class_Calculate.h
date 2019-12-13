//
// Created by Bowen Zhao on 9/26/19.
//
#pragma once
#include "Class_SpecFunc.h"
#include "Config.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace AC{
    class Class_Calculate {
    private:
        // Some constants will used in the calculation

        // Store the spectral function information
        Class_SpecFunc SpecFunc;

        // Store the mappng: grid_omega -> real omega values
        // argument: index_omega
        std::vector<AC::Type_ValReal> Grid_Omega;
        // Store the mapping: grid_tau -> real tau values
        // argument: index_tau
        std::vector<AC::Type_ValReal> Grid_Tau;
        // Store the mapping: grid_tau -> real G values
        // argument: index_tau
        // Store the mapping: grid_tau, grid_tau -> cov G values
        std::vector<AC::Type_ValReal> Array_G;
        std::vector<std::vector<AC::Type_ValReal>> Array_CovG;
        // Store the pre-cal kernel function:
        // argument: (index:\tau, index:\omega)
        std::vector<std::vector<AC::Type_ValReal>> Mat_Kernel;

        // Model Parameters
        Type_ValReal Val_Beta;

    public:
        Class_Calculate() = default;
        ~Class_Calculate() = default;
        explicit Class_Calculate(Type_ValReal _Min_Omega, Type_ValReal _Max_Omega, int _Num_DivideOmega,
                                 int _Num_DeltaFunc,
                                 std::string _file_G, std::string _file_Cov, Type_ValReal _Val_Beta);

        Type_ValReal Func_Kernel(Type_ValReal _Val_Tau, Type_ValReal _Val_Omega);

    };
}

