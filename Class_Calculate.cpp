//
// Created by Bowen Zhao on 9/26/19.
//

#include "Class_Calculate.h"

AC::Class_Calculate::Class_Calculate(Type_ValReal _min_Omega, Type_ValReal _max_Omega, int _num_DivideOmega,
                                     int _num_DeltaFunc,
                                     const std::string &_file_G, const std::string &_file_Cov,
                                     const std::string &_file_Output_para, const std::string &_file_Output,
                                     int _num_SpectralBins,
                                     Type_ValReal _val_Beta, Type_ValReal _val_WeightLeading)
        :
        specFunc(_num_DeltaFunc, _num_DivideOmega, _val_WeightLeading),
        specFunc_noLarge(_num_DeltaFunc, _num_DivideOmega, 0.),
        Grid_Omega(_num_DivideOmega + 1),
        Val_Beta(_val_Beta),
        Val_Theta(START_THETA),
        Num_SpectralBins(_num_SpectralBins),
        Str_Output_Para(_file_Output_para),
        Str_Output(_file_Output),
        Val_WeightLeading(_val_WeightLeading)
        {
    // Initialize grid_Omega
    Type_ValReal GridLength_Omega = (_max_Omega - _min_Omega) / _num_DivideOmega;

    for (int index_Omega = 0; index_Omega != _num_DivideOmega + 1; ++index_Omega) {
        Grid_Omega(index_Omega) = (_min_Omega + index_Omega * GridLength_Omega);
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
    this->specFunc.Cal_G_tilde(array_G_tilde, Mat_Kernel, Array_G);
    this->specFunc_noLarge.Cal_G_tilde(array_G_tilde_noLarge, Mat_Kernel, Array_G);

    // Current chi2 and current minimum chi2
    Val_chi2 = Cal_chi2(array_G_tilde_noLarge);
    Val_min_chi2 = Val_chi2;
}

AC::Type_ValReal AC::Class_Calculate::Func_Kernel(Type_ValReal _val_Tau, Type_ValReal _val_Omega){
    return std::exp(- _val_Tau * _val_Omega);
//    return 1./Val_PI * (std::exp(- _val_Tau * _val_Omega) + std::exp(- (Val_Beta - _val_Tau) * _val_Omega)) / (1 + std::exp(- Val_Beta * _val_Omega));
}

AC::Type_ValReal AC::Class_Calculate::Cal_chi2(const arma::Col<AC::Type_ValReal> &_array_G_tilde) {
    Type_ValReal temp_chi2 = arma::accu(_array_G_tilde % _array_G_tilde % Array_sig2inv);

    return temp_chi2;
}



void AC::Class_Calculate::Equilibrium(int _num_Steps,
                                      int _num_Bins,
                                      Class_SpecFunc & _SpecFunc,
                                      arma::Col<AC::Type_ValReal> & _array_G_tilde) {
    for (int index_Bin = 0; index_Bin != _num_Bins; ++index_Bin) {
        int num_Accept_Small = 0;
        int num_Accept_Large = 0;
        for (int index_Step = 0; index_Step != _num_Steps; ++index_Step) {
            if (Update_SmallDeltaFunc(1, _SpecFunc, _array_G_tilde)) {
                ++num_Accept_Small;
            }
            if (Update_LargeDeltaFunc(_SpecFunc, _array_G_tilde)) {
                ++num_Accept_Large;
            }
        }

        // Change the omega steps
        double ratio_Accept_Small = static_cast<double>(num_Accept_Small) / _num_Steps;
        double ratio_Accept_Large = static_cast<double>(num_Accept_Large) / _num_Steps;

        // Change small
        if (ratio_Accept_Small >= UPPER_ACCEPT) {
            _SpecFunc.Change_SmallStep(RATIO_STEPCHANGE);
        }
        else if (ratio_Accept_Small < LOWER_ACCEPT){
            _SpecFunc.Change_SmallStep(1. / RATIO_STEPCHANGE);
        }
        // Change Large
        if (ratio_Accept_Large >= UPPER_ACCEPT) {
            _SpecFunc.Change_LargeStep(RATIO_STEPCHANGE);
        }
        else if (ratio_Accept_Large < LOWER_ACCEPT){
            _SpecFunc.Change_LargeStep(1. / RATIO_STEPCHANGE);
        }
    }
}


void AC::Class_Calculate::Equilibrium_noLarge(int _num_Steps,
                                              int _num_Bins,
                                              Class_SpecFunc & _SpecFunc,
                                              arma::Col<AC::Type_ValReal> & _array_G_tilde) {
    for (int index_Bin = 0; index_Bin != _num_Bins; ++index_Bin) {
        int num_Accept_Small = 0;
        for (int index_Step = 0; index_Step != _num_Steps; ++index_Step) {
            if (Update_SmallDeltaFunc(1, _SpecFunc, _array_G_tilde)) {
                ++num_Accept_Small;
            }
        }

        // Change the omega steps
        double ratio_Accept_Small = static_cast<double>(num_Accept_Small) / _num_Steps;

        // Change small
        if (ratio_Accept_Small >= UPPER_ACCEPT) {
            _SpecFunc.Change_SmallStep(RATIO_STEPCHANGE);
        }
        else if (ratio_Accept_Small < LOWER_ACCEPT){
            _SpecFunc.Change_SmallStep(1. / RATIO_STEPCHANGE);
        }
    }
}


void AC::Class_Calculate::Anneal_Theta() {
    // Annealing
    std::vector<AC::Type_AnnealData> Array_AnnealData;
    for (int index_Anneal = 0; index_Anneal != MAX_ANNEALING; ++index_Anneal) {
        Equilibrium_noLarge(EQ_STEPS, EQ_BINS, specFunc_noLarge, array_G_tilde_noLarge);
        Type_ValReal current_chi2 = 0.;
        for (int index_Step = 0; index_Step != NUM_MEASURE_ANNEALING; ++index_Step) {
            Update_SmallDeltaFunc(1, specFunc_noLarge, array_G_tilde_noLarge);
            Update_SmallDeltaFunc(2, specFunc_noLarge, array_G_tilde_noLarge);
            current_chi2 += Val_chi2;
        }
        current_chi2 /= NUM_MEASURE_ANNEALING;

        Array_AnnealData.emplace_back(current_chi2, Val_Theta);

        if ((current_chi2 - Val_min_chi2) < THRESHOD_CHI2){
            break;
        }

        std::cout << index_Anneal << "\t" << Val_Theta <<"\t" << current_chi2 / Array_G.size()<< std::endl;

        Val_Theta /= RATIO_TEMPCHANGE;
    }

    // Get it back
    for (auto ptr_AnnealData = Array_AnnealData.cend(); ptr_AnnealData != Array_AnnealData.cbegin(); --ptr_AnnealData) {
        auto temp_chi2 = std::get<0>(*ptr_AnnealData);
//        if (temp_chi2 > Val_min_chi2 + std::sqrt(2 * Val_min_chi2)){
        if (temp_chi2 > 2 * Val_min_chi2){
            Val_Theta = std::get<1>(*ptr_AnnealData);

            std::cout << Val_Theta <<"\t" << temp_chi2 / Array_G.size()<< std::endl;
            break;
        }
    }
}

void AC::Class_Calculate::Measure_Spectral() {
    Type_Spectral temp_Spectral(Num_SpectralBins, arma::fill::zeros);
    specFunc.Measure_Spectral(temp_Spectral, Grid_Omega.size() - 1);

    Array_Spectral.AppendValue(temp_Spectral);
    Mea_Chi2.AppendValue(Val_chi2);
}

void AC::Class_Calculate::Write_Para() {
    auto sum_chi2 = Mea_Chi2.Get_AveValue();
    Mea_Chi2.ClearValue();

    std::ofstream outfile;
    outfile.open(Str_Output_Para, std::ofstream::out);
    outfile << Val_WeightLeading << "\t" << Val_Theta <<"\t" << sum_chi2 / Array_G.size()<< std::endl;

    outfile.close();

}

void AC::Class_Calculate::WriteBin_Spectral(int _index_bin) {
    auto ave_Array_Spectral = Array_Spectral.Get_AveValue();

    std::ofstream outfile;
    if (_index_bin == 0) {
        outfile.open(Str_Output, std::ofstream::out);
    }
    else {
        outfile.open(Str_Output, std::ofstream::out | std::ofstream::app);
    }

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

    std::cout << _index_bin << "\t" << Val_chi2 / Array_G.size() << std::endl;

}

void AC::Class_Calculate::CleanBin_Spectral() {
    std::ofstream outfile;
    outfile.open(Str_Output, std::ofstream::out);
    outfile.close();
}

bool AC::Class_Calculate::Update_SmallDeltaFunc(
        int _num_UpdateDeltaFunc,
        Class_SpecFunc & _SpecFunc,
        arma::Col<AC::Type_ValReal> & _array_G_tilde
        ) {
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
    static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(0., 1.);

    auto temp_Array_G_tilde = _array_G_tilde;
    std::vector<AC::Type_UpdateInfo_Small> temp_array_info;

    // Update small DeltaFunc
    for (int index_DeltaFunc = 0; index_DeltaFunc != _num_UpdateDeltaFunc + 1; ++index_DeltaFunc) {
        auto tuple_PosAmp = _SpecFunc.Get_RandomSmallDelta();

        auto pos_old = std::get<1>(tuple_PosAmp);
        auto pos_new = std::get<2>(tuple_PosAmp);
        auto amp = _SpecFunc.Get_Val_SmallAmp();

        const auto pos_LargeDelta = _SpecFunc.Get_Pos_LargeDelta();

        if ((pos_new < pos_LargeDelta) || (pos_new >= Grid_Omega.size())){
            return false;
        }

        temp_Array_G_tilde += ( Mat_Kernel.col(pos_new) - Mat_Kernel.col(pos_old)) * amp;
        temp_array_info.emplace_back(tuple_PosAmp);
    }


    // Decide if we need to update
    auto temp_chi2 = Cal_chi2(temp_Array_G_tilde);

    if (temp_chi2 < Val_chi2) {
        _SpecFunc.Update_SmallDeltaFunc(temp_array_info);
        _array_G_tilde = temp_Array_G_tilde;
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
            _SpecFunc.Update_SmallDeltaFunc(temp_array_info);
            _array_G_tilde = temp_Array_G_tilde;
            Val_chi2 = temp_chi2;

            return true;
        }
        else{
            return false;
        }
    }
}

bool AC::Class_Calculate::Update_LargeDeltaFunc(
        Class_SpecFunc & _SpecFunc,
        arma::Col<AC::Type_ValReal> & _array_G_tilde
        ) {
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
    static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(0., 1.);

    auto temp_Array_G_tilde = _array_G_tilde;
    std::vector<AC::Type_UpdateInfo_Small> temp_array_info;

    // Update Large DeltaFunc
    auto tuple_PosAmp_Large = _SpecFunc.Get_RandomLargePos();

    auto pos_old = std::get<0>(tuple_PosAmp_Large);
    auto pos_new = std::get<1>(tuple_PosAmp_Large);
    auto amp_Large = _SpecFunc.Get_Val_LargeAmp();

    const auto LeastPos_SmallDelta = _SpecFunc.Get_LeastPos_SmallDelta();
    if ((pos_new >= Grid_Omega.size()) || (pos_new < 0)) {
        return false;
    }


    // Change of lower boundary
    if (pos_new > LeastPos_SmallDelta) {

        _SpecFunc.Check_Bound(pos_new, temp_array_info, temp_Array_G_tilde);
        for (const auto &ele: temp_array_info) {
            temp_Array_G_tilde += ( Mat_Kernel.col(std::get<2>(ele)) - Mat_Kernel.col(std::get<1>(ele))) * _SpecFunc.Get_Val_SmallAmp();
        }

    }
    temp_Array_G_tilde += ( Mat_Kernel.col(pos_new) - Mat_Kernel.col(pos_old)) * amp_Large;

    // Decide if we need to update
    auto temp_chi2 = Cal_chi2(temp_Array_G_tilde);

    if (temp_chi2 < Val_chi2) {
        if (!temp_array_info.empty()) {
            _SpecFunc.Update_SmallDeltaFunc(temp_array_info);
            _SpecFunc.Set_LeastPos_SmallDelta(pos_new);
        }
        _SpecFunc.Update_largeDeltaFunc(tuple_PosAmp_Large);
        _array_G_tilde = temp_Array_G_tilde;
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
            if (!temp_array_info.empty()) {
                _SpecFunc.Update_SmallDeltaFunc(temp_array_info);
                _SpecFunc.Set_LeastPos_SmallDelta(pos_new);
            }
            _SpecFunc.Update_largeDeltaFunc(tuple_PosAmp_Large);
            _array_G_tilde = temp_Array_G_tilde;
            Val_chi2 = temp_chi2;

            return true;
        }
        else{
            return false;
        }
    }
}

void AC::Class_Calculate::Setup_Measure() {
    Val_chi2 = Cal_chi2(array_G_tilde);
    auto temp_Theta = Val_Theta;
    Val_Theta = START_THETA;

    while (temp_Theta < Val_Theta) {
        Equilibrium(EQ_STEPS, EQ_BINS, specFunc, array_G_tilde);
        for (int index_Step = 0; index_Step != NUM_MEASURE_ANNEALING; ++index_Step) {
            Update_SmallDeltaFunc(1, specFunc, array_G_tilde);
            Update_SmallDeltaFunc(2, specFunc, array_G_tilde);
            Update_LargeDeltaFunc(specFunc, array_G_tilde);
        }
        Val_Theta /= RATIO_TEMPCHANGE;
        std::cout << Val_Theta <<"\t" << Val_chi2 / Array_G.size()<< std::endl;

    }

}
