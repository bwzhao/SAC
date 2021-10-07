//
// Created by Bowen Zhao on 9/26/19.
//

#include "Class_SpecFunc.h"

AC::Class_SpecFunc::Class_SpecFunc(int _num_SmallDeltaFunc, int _num_DivideOmega, Type_ValReal _weight_LeadingOmega) :
        Num_SmallDeltaFunc(_num_SmallDeltaFunc),
        Array_Pos_SmallDelta(_num_SmallDeltaFunc),
        Val_SmallAmp(G_0 * (1 - _weight_LeadingOmega) / (_num_SmallDeltaFunc)),
        Val_LargeAmp(G_0 * _weight_LeadingOmega),
        Val_SmallStep(static_cast<int>(static_cast<double>(_num_DivideOmega + 1) / (_num_SmallDeltaFunc - 1) * 10)),
        Val_LargeStep(Val_SmallStep)
{
    // Initialize the A(w) function
    // Use a uniform distribution
    for (int index_DeltaFunc = 0; index_DeltaFunc != _num_SmallDeltaFunc; ++index_DeltaFunc) {
        auto temp_pos = static_cast<int>(static_cast<double>(_num_DivideOmega + 1) / (_num_SmallDeltaFunc - 1)) * index_DeltaFunc;
        Array_Pos_SmallDelta(index_DeltaFunc) = temp_pos;
    }
    LeastPos_SmallDelta = 0;
    Pos_LargeDelta = 0;
}

void AC::Class_SpecFunc::Cal_G_tilde(arma::Col<AC::Type_ValReal> &_array_GTilde,
                                     const arma::Mat<AC::Type_ValReal> &_mat_Kernel,
                                     const arma::Col<AC::Type_ValReal> &_array_G) {
    _array_GTilde = -_array_G;
    for (int index_Tau = 0; index_Tau != _array_GTilde.size(); ++index_Tau) {
        for (int index_DeltaFunc = 0; index_DeltaFunc != Num_SmallDeltaFunc; ++index_DeltaFunc) {
            _array_GTilde(index_Tau) += _mat_Kernel(index_Tau, Array_Pos_SmallDelta(index_DeltaFunc)) * Val_SmallAmp;
        }
        _array_GTilde(index_Tau) += _mat_Kernel(index_Tau, Pos_LargeDelta) * Val_LargeAmp;
    }
}


//AC::Type_UpdateInfo_Small AC::Class_SpecFunc::Get_RandomSmallDelta() {
//    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
//    static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(-1, 1.);
//    static std::uniform_int_distribution<int> Ran_Int(0., Num_SmallDeltaFunc - 1);
//
//    int which_index = Ran_Int(engine);
//
//    int pos_old = Array_Pos_SmallDelta(which_index);
//    int pos_new = Array_Pos_SmallDelta(which_index) + int(Ran_Double(engine) * Val_SmallStep);
//
//    double amp = Array_Amp(which_index);
//
//    return std::make_tuple(which_index, pos_old, pos_new, amp);
//}


void AC::Class_SpecFunc::Measure_Spectral(Type_Spectral &_array_Spectral, int _len_GridOmega) {
    int Binlength = _len_GridOmega / _array_Spectral.size();

    // Small Delta-func
    for (int index_DeltaFunc = 0; index_DeltaFunc != Num_SmallDeltaFunc; ++index_DeltaFunc) {
        auto index_bin = Array_Pos_SmallDelta[index_DeltaFunc] / Binlength;

        if (index_bin != _array_Spectral.size()){
            _array_Spectral[index_bin] += Val_SmallAmp;
        }
        // Right boundary belongs to the last bin
        else{
            _array_Spectral[index_bin - 1] += Val_SmallAmp;
        }
    }

    // Large Delta-func
    auto index_bin = Pos_LargeDelta / Binlength;

    if (index_bin != _array_Spectral.size()){
        _array_Spectral[index_bin] += Val_LargeAmp;
    }
        // Right boundary belongs to the last bin
    else{
        _array_Spectral[index_bin - 1] += Val_LargeAmp;
    }
}

