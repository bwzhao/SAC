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
#include <armadillo>

namespace AC{
    class Class_Calculate {
    private:
        // Some constants will used in the calculation

        // Store the spectral function information
        Class_SpecFunc SpecFunc;

        // Store the mappng: grid_omega -> real omega values
        // argument: index_omega
        arma::Col<AC::Type_ValReal> Grid_Omega;

        // Store the mapping: grid_tau -> real tau values
        // argument: index_tau
        arma::Col<AC::Type_ValReal> Grid_Tau;

        // Store the mapping: grid_tau -> real G values
        // argument: index_tau
        // Store the mapping: grid_tau, grid_tau -> cov G values
        arma::Col<AC::Type_ValReal> Array_G;
        arma::Mat<AC::Type_ValReal> Mat_CovG;
        // Inverse sigma
        arma::Col<AC::Type_ValReal> Array_sig2inv;

        // Store the pre-cal kernel function:
        // argument: (index:\tau, index:\omega)
        arma::Mat<AC::Type_ValReal> Mat_Kernel;

        // Store the differencce between estimation of G from A(w) and the real value of G
        // G_Tilde = estimation - G
        // argument: (index:\tau)
        arma::Col<AC::Type_ValReal> Array_G_tilde;
        // Current value of chi2
        Type_ValReal Val_chi2;

        // Model Parameters
        Type_ValReal Val_Beta;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Parameters using the the calculation
        // Sampling temperature
        Type_ValReal Val_Theta;



    public:
        Class_Calculate() = default;
        ~Class_Calculate() = default;
        explicit Class_Calculate(Type_ValReal _Min_Omega, Type_ValReal _Max_Omega, int _Num_DivideOmega,
                                 int _Num_DeltaFunc,
                                 std::string _file_G, std::string _file_Cov, Type_ValReal _Val_Beta);

        Type_ValReal Func_Kernel(Type_ValReal _Val_Tau, Type_ValReal _Val_Omega);

        AC::Type_ValReal Cal_chi2(const arma::Col<AC::Type_ValReal> &_Array_G_tilde);

        void Update_One();
    };
}

