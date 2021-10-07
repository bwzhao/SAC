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
#include "Class_Spectral_Measurement.h"

namespace AC{

    class Class_Calculate {

        // Some constants will used in the calculation
    public:
        // Store the spectral function information
        Class_SpecFunc specFunc;
        Class_SpecFunc specFunc_noLarge;

        // Store the differencce between estimation of G from A(w) and the real value of G
        // G_Tilde = estimation - G
        // argument: (index:\tau)
        arma::Col<AC::Type_ValReal> array_G_tilde;
        arma::Col<AC::Type_ValReal> array_G_tilde_noLarge;

        // Parameters using the the calculation
        // Sampling temperature
        Type_ValReal Val_Theta;
    private:
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

        // Current value of chi2
        Type_ValReal Val_chi2;
        Type_ValReal Sum_chi2;
        Type_ValReal Val_min_chi2;

        // Model Parameters
        Type_ValReal Val_Beta;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        Type_ValReal Val_WeightLeading;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Parameters for Measurement
        int Num_SpectralBins;
        Class_Spectral_Measurement<AC::Type_ValReal> Array_Spectral;

        std::string Str_Output_Para;
        std::string Str_Output;

        Class_Data_Measurement<AC::Type_ValReal> Mea_Chi2;

    public:
        Class_Calculate() = default;
        ~Class_Calculate() = default;
        explicit Class_Calculate(Type_ValReal _min_Omega, Type_ValReal _max_Omega, int _num_DivideOmega,
                                 int _num_DeltaFunc,
                                 const std::string &_file_G, const std::string &_file_Cov,
                                 const std::string &_file_Output_para, const std::string &_file_Output,
                                 int _num_SpectralBins,
                                 Type_ValReal _val_Beta,
                                 Type_ValReal _val_WeightLeading);

        // Kernel Function
        Type_ValReal Func_Kernel(Type_ValReal _val_Tau, Type_ValReal _val_Omega);

        // Calculate chi2 based on current G_tilde
        Type_ValReal Cal_chi2(const arma::Col<AC::Type_ValReal> &_array_G_tilde);

        // Update process
        // Just move one delta function at a time
        bool Update_SmallDeltaFunc(int _num_UpdateDeltaFunc,
                                   Class_SpecFunc & _SpecFunc,
                                   arma::Col<AC::Type_ValReal> & _array_G_tilde);
        bool Update_LargeDeltaFunc(Class_SpecFunc & _SpecFunc,
                                   arma::Col<AC::Type_ValReal> & _array_G_tilde);

        // Equilibrium
        void Equilibrium(int _num_Steps, int _num_Bins,
                         Class_SpecFunc &_SpecFunc,
                         arma::Col<AC::Type_ValReal> & _array_G_tilde);
        void Equilibrium_noLarge(int _num_Steps, int _num_Bins,
                                 Class_SpecFunc &_SpecFunc,
                                 arma::Col<AC::Type_ValReal> & _array_G_tilde);

        // Annealing theta to get the optimal value
        void Anneal_Theta();

        // Measure specture
        void Measure_Spectral();
        void WriteBin_Spectral(int _index_bin);

        void Measure_Para();
        void Write_Para();

        void CleanBin_Spectral();
        void Setup_Measure();
    };
}

