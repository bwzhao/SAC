//
// Created by Bowen Zhao on 9/26/19.
//

#include "Class_Calculate.h"

AC::Class_Calculate::Class_Calculate(Type_ValReal _Min_Omega, Type_ValReal _Max_Omega, int _Num_OmegaGrid,
                                     int _Num_DeltaFunc, std::string _file_G, std::string _file_Cov,
                                     std::string _file_Output, int _Num_SpectralBins, Type_ValReal _Val_Beta)
        :
        SpecFunc(_Num_DeltaFunc, _Num_OmegaGrid),
        Val_Beta(_Val_Beta),
        Val_Theta(START_THETA),
        Grid_Omega(_Num_OmegaGrid + 1),
        Num_SpectralBins(_Num_SpectralBins),
        Str_Output(_file_Output)
        {
    // Initialize grid_Omega
    Type_ValReal GridLength_Omega = (_Max_Omega - _Min_Omega) / _Num_OmegaGrid;
    for (int index_Omega = 0; index_Omega != _Num_OmegaGrid + 1; ++index_Omega) {
        Grid_Omega(index_Omega) = (_Min_Omega + index_Omega * GridLength_Omega);
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
        Type_ValReal temp_Gerr = 0.;

        file_G >> temp_tau;
        if (file_G.peek() == EOF){
            break;
        }
        file_G >> temp_G;
        if (file_G.peek() == EOF){
            break;
        }
        file_G >> temp_Gerr;
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

    arma::eig_sym(eigval, eigvec, Mat_CovG.i());

    Array_G = eigvec.t() * Array_G;
    Mat_Kernel =  eigvec.t() * Mat_Kernel;
    Array_sig2inv = eigval;

    // Initialize G_tilde
    this->SpecFunc.Cal_G_tilde(Array_G_tilde, Mat_Kernel, Array_G);

    // Current chi2 and current minimum chi2
    Val_chi2 = Cal_chi2(Array_G_tilde);
    Val_min_chi2 = Val_chi2;
}

AC::Type_ValReal AC::Class_Calculate::Func_Kernel(Type_ValReal _Val_Tau, Type_ValReal _Val_Omega){
    return std::exp(- _Val_Tau * _Val_Omega);
//    return 1./Val_PI * (std::exp(- _Val_Tau * _Val_Omega) + std::exp(- (Val_Beta - _Val_Tau) * _Val_Omega)) / (1 + std::exp(- Val_Beta * _Val_Omega));
}

AC::Type_ValReal AC::Class_Calculate::Cal_chi2(const arma::Col<AC::Type_ValReal> &_Array_G_tilde) {
    Type_ValReal temp_chi2 = arma::accu(_Array_G_tilde % _Array_G_tilde % Array_sig2inv);

    return temp_chi2;
}

bool AC::Class_Calculate::Update_DeltaFunc(int _Num_DeltaFunc) {
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
    static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(0., 1.);

    auto temp_Array_G_tilde = Array_G_tilde;
    std::vector<AC::Type_UpdateInfo> temp_array_info;

    for (int index_DeltaFunc = 0; index_DeltaFunc != _Num_DeltaFunc; ++index_DeltaFunc) {
        auto tuple_PosAmp = SpecFunc.Get_RandomPos();
        auto index = std::get<0>(tuple_PosAmp);
        auto pos_old = std::get<1>(tuple_PosAmp);
        auto pos_new = std::get<2>(tuple_PosAmp);
        auto amp = std::get<3>(tuple_PosAmp);

        if ((pos_new < 0) || (pos_new >= Grid_Omega.size())){
            return false;
        }

        temp_Array_G_tilde += ( Mat_Kernel.col(pos_new) - Mat_Kernel.col(pos_old)) * amp;
        temp_array_info.emplace_back(tuple_PosAmp);
    }


    // Decide if we need to update
    auto temp_chi2 = Cal_chi2(temp_Array_G_tilde);

    if (temp_chi2 < Val_chi2) {
        SpecFunc.UpdateOne(temp_array_info);
        Array_G_tilde = temp_Array_G_tilde;
        Val_chi2 = temp_chi2;

        // Update minimnm chi2
        if (temp_chi2 < Val_min_chi2) {
            Val_min_chi2 = temp_chi2;
        }

        return true;
    }
    else{
        auto ratio = std::exp(-(temp_chi2 - Val_chi2) / (2. * Val_Theta));
        // If accepting the change
        if (ratio > Ran_Double(engine)) {
            SpecFunc.UpdateOne(temp_array_info);
            Array_G_tilde = temp_Array_G_tilde;
            Val_chi2 = temp_chi2;

            return true;
        }
        else{
            return false;
        }
    }
}

void AC::Class_Calculate::Equilibrium(int _Num_Steps, int _Num_Bins) {
    for (int index_Bin = 0; index_Bin != _Num_Bins; ++index_Bin) {
        int num_Accept = 0;
        for (int index_Step = 0; index_Step != _Num_Steps; ++index_Step) {
            if (Update_DeltaFunc(1)) {
                ++num_Accept;
            }
        }
        // Change the omega steps
        double ratio_Accept = static_cast<double>(num_Accept) / _Num_Steps;
        if (ratio_Accept >= UPPER_ACCEPT) {
            SpecFunc.Change_OmegaStep(RATIO_STEPCHANGE);
        }
        else if (ratio_Accept < LOWER_ACCEPT){
            SpecFunc.Change_OmegaStep(1. / RATIO_STEPCHANGE);
        }
    }
}

AC::Type_ValReal AC::Class_Calculate::Measure_chi2(int _Num_Steps) {
    Type_ValReal temp_val_chi2 = 0.;
    for (int index_Step = 0; index_Step != _Num_Steps; ++index_Step) {
        Update_DeltaFunc(1);
        Update_DeltaFunc(2);
        temp_val_chi2 += Val_chi2;
    }
    temp_val_chi2 /= _Num_Steps;

    return temp_val_chi2;
}

void AC::Class_Calculate::Anneal_Theta() {
    // Annealing
    std::vector<AC::Type_AnnealData> Array_AnnealData;
    for (int index_Anneal = 0; index_Anneal != MAX_ANNEALING; ++index_Anneal) {
        Equilibrium(EQ_STEPS, EQ_BINS);
        Type_ValReal current_chi2 = Measure_chi2(NUM_MEASURE_ANNEALING);
        Array_AnnealData.emplace_back(current_chi2, Val_Theta, SpecFunc);

        if ((current_chi2 - Val_min_chi2) < THRESHOD_CHI2){
            break;
        }
        std::cout << index_Anneal << "\t" << Val_Theta <<"\t" << current_chi2 / Array_G.size()<< std::endl;

        Val_Theta /= RATIO_TEMPCHANGE;
    }

    // Get it back
    for (auto ptr_AnnealData = Array_AnnealData.cend(); ptr_AnnealData != Array_AnnealData.cbegin(); --ptr_AnnealData) {
        auto temp_chi2 = std::get<0>(*ptr_AnnealData);
        if (temp_chi2 > Val_min_chi2 + std::sqrt(2 * Val_min_chi2)){
            Val_Theta = std::get<1>(*ptr_AnnealData);
            SpecFunc = std::get<2>(*ptr_AnnealData);

            std::cout <<Val_Theta <<"\t" << temp_chi2 / Array_G.size()<< std::endl;
            break;
        }
    }
}

void AC::Class_Calculate::Measure_Spectral() {
    Type_Spectral temp_Spectral(Num_SpectralBins, arma::fill::zeros);
    SpecFunc.Measure_Spectral(temp_Spectral, Grid_Omega.size() - 1);

    Array_Spectral.AppendValue(temp_Spectral);
}

void AC::Class_Calculate::WriteBin_Spectral() {
    auto ave_Array_Spectral = Array_Spectral.Get_AveValue();

    std::ofstream outfile;
    outfile.open(Str_Output, std::ofstream::out | std::ofstream::app);
//    outfile.setf(std::ios::fixed);
//    outfile.precision(12);

    for (int index_Spectral = 0; index_Spectral != ave_Array_Spectral.size() - 1; ++index_Spectral) {
        auto Omega_midpoint = 0.5 * (Grid_Omega[index_Spectral * Grid_Omega.size() / Num_SpectralBins] + Grid_Omega[static_cast<int>((index_Spectral + 1) * Grid_Omega.size() / Num_SpectralBins)]);
        auto Length_Bin = (Grid_Omega[static_cast<int>((index_Spectral + 1) * Grid_Omega.size() / Num_SpectralBins)] - Grid_Omega[static_cast<int>((index_Spectral) * Grid_Omega.size() / Num_SpectralBins)]);

        outfile << Omega_midpoint << "\t"
                << ave_Array_Spectral[index_Spectral] / Length_Bin << std::endl;
    }
    outfile.close();

    Array_Spectral.ClearValue();
}

void AC::Class_Calculate::CleanBin_Spectral() {
    std::ofstream outfile;
    outfile.open(Str_Output, std::ofstream::out);
    outfile.close();
}
