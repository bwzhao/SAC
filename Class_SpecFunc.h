//
// Created by Bowen Zhao on 9/26/19.
//
#pragma once

#include <vector>
#include "Config.h"
#include <armadillo>
#include <random>

namespace AC{
    class Class_SpecFunc {
    private:
        // Number of the total delta function
        int Num_DeltaFunc;

        // Store the position of the delta function
        arma::Col<int> Array_Pos;

        // Store the amplitude of the delta function
        arma::Col<Type_ValReal> Array_Amp;

        // Update step for Omega
        Type_ValReal Val_OmegaStep;

    public:
        Class_SpecFunc() = default;
        explicit Class_SpecFunc(int _Num_DeltaFunc, int _Num_DivideOmega);

        int size() const{
            return Num_DeltaFunc;
        }

        void Cal_G_tilde(arma::Col<AC::Type_ValReal> &_Array_GTilde, const arma::Mat<AC::Type_ValReal> &_Mat_Kernel,
                         const arma::Col<AC::Type_ValReal> &_Array_G);

        void UpdateOne(const std::tuple<int, int> & _tuple_Data);

        std::tuple<int, int, int, AC::Type_ValReal> Get_RandomPos();

        void Change_OmegaStep(Type_ValReal _ratio){
            Val_OmegaStep *= _ratio;
        }

        int Get_Num_DeltaFunc(){
            return Num_DeltaFunc;
        }

        void Measure_Spectral(Type_Spectral &_array_Spectral, int _Num_DivideOmega);

    };
}



