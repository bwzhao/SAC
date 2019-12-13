//
// Created by Bowen Zhao on 9/26/19.
//

#include "Class_Calculate.h"

AC::Class_Calculate::Class_Calculate(Type_ValReal _Min_Omega, Type_ValReal _Max_Omega, int _Num_DivideOmega,
                                     int _Num_DeltaFunc,
                                     std::string _file_G, std::string _file_Cov, Type_ValReal _Val_Beta):
        Val_Beta(_Val_Beta){
    // Initialize grid_Omega
    Type_ValReal Delta_Omega = (_Max_Omega - _Min_Omega) / _Num_DivideOmega;
    for (int index_Omega = 0; index_Omega != _Num_DivideOmega + 1; ++index_Omega) {
        Grid_Omega.emplace_back(_Min_Omega + index_Omega * Delta_Omega);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize val_G and Grid_Tau
    std::ifstream file_G(_file_G);
    if(!file_G.good()){
        std::cout << "There is no data file!" << std::endl;
        file_G.close();
        exit(-1);
    }
    while (file_G.peek() != EOF) {
        Type_ValReal temp_tau = 0.;
        Type_ValReal temp_G = 0.;

        file_G >> temp_tau;
        if (file_G.peek() == EOF){
            break;
        }
        file_G >> temp_G;
        if (file_G.peek() == EOF){
            break;
        }

        Grid_Tau.emplace_back(temp_tau);
        Array_G.emplace_back(temp_G);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize Array_CovG
    std::ifstream file_CovG(_file_Cov);
    if(!file_CovG.good()){
        std::cout << "There is no data file!" << std::endl;
        file_CovG.close();
        exit(-1);
    }
    for (int index_row = 0; index_row != Grid_Tau.size(); ++index_row) {
        std::vector<AC::Type_ValReal> temp_row;
        for (int index_column = 0; index_column != Grid_Tau.size(); ++index_column) {
            Type_ValReal temp_CovEle = 0;
            file_CovG >> temp_CovEle;
            temp_row.emplace_back(temp_CovEle);
        }
        Array_CovG.emplace_back(temp_row);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize Kernel
    for (int index_Tau = 0; index_Tau != Grid_Tau.size(); ++index_Tau) {
        std::vector<AC::Type_ValReal> temp_row;
        for (int index_Omega = 0; index_Omega != Grid_Omega.size(); ++index_Omega) {
            temp_row.emplace_back(Func_Kernel(Grid_Tau[index_Tau], Grid_Omega[index_Omega]));
        }
        Mat_Kernel.emplace_back(temp_row);
    }
}

AC::Type_ValReal AC::Class_Calculate::Func_Kernel(Type_ValReal _Val_Tau, Type_ValReal _Val_Omega){
    return 1./Val_PI * std::exp(- _Val_Tau * _Val_Omega) / (1 + std::exp(Val_Beta * _Val_Tau));
}
