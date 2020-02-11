//
// Created by Bowen Zhao on 9/26/19.
//

#include "Class_Calculate.h"

AC::Class_Calculate::Class_Calculate(Type_ValReal _Min_Omega, Type_ValReal _Max_Omega, int _Num_DivideOmega,
                                     int _Num_DeltaFunc,
                                     std::string _file_G, std::string _file_Cov, Type_ValReal _Val_Beta):
        SpecFunc(_Num_DeltaFunc, _Num_DivideOmega),
        Val_Beta(_Val_Beta),
        Val_Theta(10.),
        Grid_Omega(_Num_DivideOmega + 1)
        {
    // Initialize grid_Omega
    Type_ValReal Delta_Omega = (_Max_Omega - _Min_Omega) / _Num_DivideOmega;
    for (int index_Omega = 0; index_Omega != _Num_DivideOmega + 1; ++index_Omega) {
        Grid_Omega(index_Omega) = (_Min_Omega + index_Omega * Delta_Omega);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize val_G and Grid_Tau
    std::ifstream file_G(_file_G);
    if(!file_G.good()){
        std::cout << "There is no data file!" << std::endl;
        file_G.close();
        exit(-1);
    }

    auto temp_Grid_Tau = std::vector<AC::Type_ValReal>();
    auto temp_Array_G = std::vector<AC::Type_ValReal>();
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

        temp_Grid_Tau.emplace_back(temp_tau);
        temp_Array_G.emplace_back(temp_G);
    }
    Array_G = decltype(Array_G)(temp_Array_G);
    Grid_Tau = decltype(Grid_Tau)(temp_Grid_Tau);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize Mat_CovG
    std::ifstream file_CovG(_file_Cov);
    if(!file_CovG.good()){
        std::cout << "There is no data file!" << std::endl;
        file_CovG.close();
        exit(-1);
    }

    Mat_CovG = decltype(Mat_CovG)(Grid_Tau.size(), Grid_Tau.size());
    for (int index_row = 0; index_row != Grid_Tau.size(); ++index_row) {
        for (int index_column = 0; index_column != Grid_Tau.size(); ++index_column) {
            Type_ValReal temp_CovEle = 0;
            file_CovG >> temp_CovEle;
            Mat_CovG(index_row, index_column) = temp_CovEle;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize Kernel
    Mat_Kernel = decltype(Mat_Kernel)(Grid_Tau.size(), Grid_Omega.size());
    for (int index_Tau = 0; index_Tau != Grid_Tau.size(); ++index_Tau) {
        for (int index_Omega = 0; index_Omega != Grid_Omega.size(); ++index_Omega) {
            Mat_Kernel(index_Tau, index_Omega) = (Func_Kernel(Grid_Tau(index_Tau), Grid_Omega(index_Omega)));
        }
    }

    // Transfer to eigen basis
    arma::Col<AC::Type_ValReal> eigval;
    arma::Mat<AC::Type_ValReal> eigvec;

    arma::eig_sym(eigval, eigvec, Mat_CovG);

    Array_G = eigvec.t() * Array_G;
    Mat_Kernel =  eigvec.t() * Mat_Kernel;
    Array_sig2inv = 1. / eigval;

    // Initialize G_tilde
    this->SpecFunc.Cal_GTilde(Array_G_tilde, Mat_Kernel, Array_G);

    Val_chi2 = Cal_chi2(Array_G_tilde);
}

AC::Type_ValReal AC::Class_Calculate::Func_Kernel(Type_ValReal _Val_Tau, Type_ValReal _Val_Omega){
    return 1./Val_PI * std::exp(- _Val_Tau * _Val_Omega) / (1 + std::exp(Val_Beta * _Val_Omega));
}

AC::Type_ValReal AC::Class_Calculate::Cal_chi2(const arma::Col<AC::Type_ValReal> &_Array_G_tilde) {
    Type_ValReal temp_chi2 = arma::accu(_Array_G_tilde % _Array_G_tilde % Array_sig2inv);

    return temp_chi2;
}

void AC::Class_Calculate::Update_One() {
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
    static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(0., 1.);

    auto tuple_PosAmp = SpecFunc.Get_RandomPos();
    auto index = std::get<0>(tuple_PosAmp);
    auto pos_old = std::get<1>(tuple_PosAmp);
    auto pos_new = std::get<2>(tuple_PosAmp);
    auto amp = std::get<3>(tuple_PosAmp);
    while (true) {
        if (pos_new > 0 && pos_new < Grid_Omega.size()){
            break;
        }
        tuple_PosAmp = SpecFunc.Get_RandomPos();
        index = std::get<0>(tuple_PosAmp);
        pos_old = std::get<1>(tuple_PosAmp);
        pos_new = std::get<2>(tuple_PosAmp);
        amp = std::get<3>(tuple_PosAmp);
    }



    auto temp_Array_G_tilde = Array_G_tilde;
    for (int index_Tau = 0; index_Tau != temp_Array_G_tilde.size(); ++index_Tau) {
        temp_Array_G_tilde(index_Tau) += (Mat_Kernel(index_Tau, pos_new) - Mat_Kernel(index_Tau, pos_old)) * amp;
    }

    auto temp_chi2 = Cal_chi2(temp_Array_G_tilde);
    auto ratio = std::exp(-(temp_chi2 - Val_chi2) / (2. * Val_Theta));
    // If accepting the change
    if (ratio > Ran_Double(engine)) {
        SpecFunc.UpdateOne(std::make_tuple(index, pos_new));
        Array_G_tilde = temp_Array_G_tilde;
        Val_chi2 = temp_chi2;
    }

}
