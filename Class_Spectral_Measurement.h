////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once

#include <map>
#include <string>
#include "Config.h"
#include <numeric>
#include <cmath>
#include <complex>
#include <armadillo>

namespace AC {
    template <typename T>
    class Class_Spectral_Measurement {
    private:
        arma::Col<T> Data;
        int Length;

    public:
        Class_Spectral_Measurement(): Length(0){};
//        Class_Spectral_Measurement(SSE::type_NumSegment _length):
//                Length(_length) {}
        ~Class_Spectral_Measurement() = default;

        // Add a measurement value to a specific quantity(key)
        void AppendValue(const arma::Col<T> &_which_value){
            if (Length == 0) {
                Data = _which_value;
            }
            else{
                Data += _which_value;
//                Data = Data + _which_value;
//                for (int index_r = 0; index_r != Data.size(); ++index_r) {
//                    Data[index_r] += _which_value[index_r];
//                }
            }
            ++Length;
        }

        arma::Col<T> Get_AveValue(){
//            std::vector<T> temp_return_vec(Data);
//            for (int index_r = 0; index_r != Data.size(); ++index_r) {
//                temp_return_vec[index_r] /= Length;
//            }
//            std::cout << Length << std::endl;
            return Data / Length;
        }

        void ClearValue(){
            Data.clear();
            Length = 0;
        }
    };
}

