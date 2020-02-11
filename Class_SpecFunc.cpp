//
// Created by Bowen Zhao on 9/26/19.
//

#include "Class_SpecFunc.h"

AC::Class_SpecFunc::Class_SpecFunc(int _Num_DeltaFunc, int _Num_DivideOmega):
        Num_DeltaFunc(_Num_DeltaFunc),
        Val_OmegaStep(static_cast<int>((_Num_DivideOmega + 1) / _Num_DeltaFunc)){
    // Initialize the A(w) function
    for (int index_DeltaFunc = 0; index_DeltaFunc != _Num_DeltaFunc; ++index_DeltaFunc) {
        Array_Pos(index_DeltaFunc) = static_cast<int>((_Num_DivideOmega + 1) / _Num_DeltaFunc * index_DeltaFunc);
        Array_Amp(index_DeltaFunc) = 1. / _Num_DeltaFunc;
    }
}

void AC::Class_SpecFunc::Cal_GTilde(arma::Col<AC::Type_ValReal> &_Array_GTilde,
                                    const arma::Mat<AC::Type_ValReal> &_Mat_Kernel,
                                    const arma::Col<AC::Type_ValReal> &_Array_G) {
    for (int index_Tau = 0; index_Tau != _Array_GTilde.size(); ++index_Tau) {
        _Array_GTilde(index_Tau) = -_Array_G(index_Tau);
        for (int index_DeltaFunc = 0; index_DeltaFunc != Num_DeltaFunc; ++Num_DeltaFunc) {
            _Array_GTilde(index_Tau)  += _Mat_Kernel(index_Tau, Array_Pos(index_DeltaFunc)) * Array_Amp(index_DeltaFunc);
        }
    }
}

void AC::Class_SpecFunc::UpdateOne(const std::tuple<int, int> & _tuple_Data) {
    const auto &index = std::get<0>(_tuple_Data);
    const auto &pos_1 = std::get<1>(_tuple_Data);

    Array_Pos(index) = pos_1;


}

std::tuple<int, int, int, AC::Type_ValReal> AC::Class_SpecFunc::Get_RandomPos() {
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
    static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(0., 1.);
    static std::uniform_int_distribution<int> Ran_Int(0., Num_DeltaFunc - 1);

    int which_index = Ran_Int(engine);

    int pos_1 = Array_Pos(which_index);
    int pos_2 = Array_Pos(which_index) + int(2 * (Ran_Double(engine) - 0.5) * Val_OmegaStep);


    double amp = Array_Amp(which_index);

    return std::make_tuple(which_index, pos_1, pos_2, amp);
}

