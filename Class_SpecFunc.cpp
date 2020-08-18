//
// Created by Bowen Zhao on 9/26/19.
//

#include "Class_SpecFunc.h"

AC::Class_SpecFunc::Class_SpecFunc(int _Num_DeltaFunc, int _GridLength_Omega):
        Num_DeltaFunc(_Num_DeltaFunc),
        Val_OmegaStep(_GridLength_Omega / _Num_DeltaFunc * 10),
        Array_Pos(_Num_DeltaFunc),
        Array_Amp(_Num_DeltaFunc)
        {
    // Initialize the A(w) function
    for (int index_DeltaFunc = 0; index_DeltaFunc != _Num_DeltaFunc; ++index_DeltaFunc) {
        auto temp_pos = static_cast<int>(static_cast<double>(_GridLength_Omega) / _Num_DeltaFunc * index_DeltaFunc);
        Array_Pos(index_DeltaFunc) = temp_pos;
        Array_Amp(index_DeltaFunc) = G_0 / _Num_DeltaFunc;
    }
}

void AC::Class_SpecFunc::Cal_G_tilde(arma::Col<AC::Type_ValReal> &_Array_GTilde,
                                     const arma::Mat<AC::Type_ValReal> &_Mat_Kernel,
                                     const arma::Col<AC::Type_ValReal> &_Array_G) {
    _Array_GTilde = -_Array_G;
    for (int index_Tau = 0; index_Tau != _Array_GTilde.size(); ++index_Tau) {
        for (int index_DeltaFunc = 0; index_DeltaFunc != Num_DeltaFunc; ++index_DeltaFunc) {
            _Array_GTilde(index_Tau)  += _Mat_Kernel(index_Tau, Array_Pos(index_DeltaFunc)) * Array_Amp(index_DeltaFunc);
        }
    }
}

void AC::Class_SpecFunc::UpdateOne(const std::vector<AC::Type_UpdateInfo> & _tuple_Data) {
    for (const auto & ele: _tuple_Data) {
        const auto &index = std::get<0>(ele);
        const auto &pos_new = std::get<2>(ele);

        Array_Pos(index) = pos_new;
    }

}

AC::Type_UpdateInfo AC::Class_SpecFunc::Get_RandomPos() {
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
    static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(-1, 1.);
    static std::uniform_int_distribution<int> Ran_Int(0., Num_DeltaFunc - 1);

    int which_index = Ran_Int(engine);

    int pos_old = Array_Pos(which_index);
    int pos_new = Array_Pos(which_index) + int(Ran_Double(engine) * Val_OmegaStep);

    double amp = Array_Amp(which_index);

    return std::make_tuple(which_index, pos_old, pos_new, amp);
}

void AC::Class_SpecFunc::Measure_Spectral(Type_Spectral &_array_Spectral, int _GridLength_Omega) {
    int Binlength = _GridLength_Omega / _array_Spectral.size();
    for (int index_DeltaFunc = 0; index_DeltaFunc != Num_DeltaFunc; ++index_DeltaFunc) {
        auto index_bin = Array_Pos[index_DeltaFunc] / Binlength;
        if (index_bin != _array_Spectral.size()){
            _array_Spectral[index_bin] += Array_Amp[index_DeltaFunc];
        }
        else{
            _array_Spectral[index_bin - 1] += Array_Amp[index_DeltaFunc];
        }

    }
}

