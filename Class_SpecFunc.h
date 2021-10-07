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
        int Num_SmallDeltaFunc;

        // Store the position of the delta function
        arma::Col<int> Array_Pos_SmallDelta;
        int LeastPos_SmallDelta;

        int Pos_LargeDelta;

        // Store the amplitude of the delta function
        Type_ValReal Val_SmallAmp;
        Type_ValReal Val_LargeAmp;

        // Update step for small delta-fcn
        Type_ValReal Val_SmallStep;
        // Update step for large delta-fcn
        Type_ValReal Val_LargeStep;

    public:
        Class_SpecFunc() = default;
        explicit Class_SpecFunc(int _num_SmallDeltaFunc, int _num_DivideOmega, Type_ValReal _weight_LeadingOmega);

        int size() const{
            return Num_SmallDeltaFunc;
        }

        int Get_Pos_LargeDelta() const{
            return Pos_LargeDelta;
        }
        int Get_LeastPos_SmallDelta() const{
            return LeastPos_SmallDelta;
        }
        int Set_LeastPos_SmallDelta(int _new_Pos){
            return LeastPos_SmallDelta = _new_Pos;
        }

        Type_ValReal Get_Val_SmallAmp() const{
            return Val_SmallAmp;
        }
        Type_ValReal Get_Val_LargeAmp() const{
            return Val_LargeAmp;
        }


        void Cal_G_tilde(arma::Col<AC::Type_ValReal> &_array_GTilde, const arma::Mat<AC::Type_ValReal> &_mat_Kernel,
                         const arma::Col<AC::Type_ValReal> &_array_G);

        void Update_SmallDeltaFunc(const std::vector<AC::Type_UpdateInfo_Small> & _array_UpdateInfo){
            for (const auto & ele: _array_UpdateInfo) {
                const auto &index = std::get<0>(ele);
                const auto &pos_new = std::get<2>(ele);

                Array_Pos_SmallDelta(index) = pos_new;

                if (pos_new < LeastPos_SmallDelta) {
                    LeastPos_SmallDelta = pos_new;
                }
            }
        }

        void Update_largeDeltaFunc(const Type_UpdateInfo_Large & _array_UpdateInfo){
            Pos_LargeDelta = std::get<1>(_array_UpdateInfo);
        }

        Type_UpdateInfo_Small Get_RandomSmallDelta(){
            static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
            static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(-1, 1.);
            static std::uniform_int_distribution<int> Ran_Int(0., Num_SmallDeltaFunc - 1);

            int which_index = Ran_Int(engine);

            int pos_old = Array_Pos_SmallDelta(which_index);
            int pos_new = Array_Pos_SmallDelta(which_index) + int(Ran_Double(engine) * Val_SmallStep);

            return std::make_tuple(which_index, pos_old, pos_new);
        }

        Type_UpdateInfo_Large Get_RandomLargePos(){
            static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr)));
            static std::uniform_real_distribution<AC::Type_ValReal> Ran_Double(-1, 1.);

            int pos_old = Pos_LargeDelta;
            int pos_new = Pos_LargeDelta + int(Ran_Double(engine) * Val_LargeStep);

            return std::make_tuple(pos_old, pos_new);
        }

        void Change_SmallStep(Type_ValReal _ratio){
            Val_SmallStep *= _ratio;
        }

        void Change_LargeStep(Type_ValReal _ratio){
            Val_LargeStep *= _ratio;
        }

        int Get_Num_DeltaFunc(){
            return Num_SmallDeltaFunc;
        }

        void Measure_Spectral(Type_Spectral &_array_Spectral, int _len_GridOmega);

        void Check_Bound(int _pos_Large, std::vector<AC::Type_UpdateInfo_Small> & _array_info, arma::Col<AC::Type_ValReal> & _array_G_tilde) const{

            for (int index_DeltaFunc = 0; index_DeltaFunc != Num_SmallDeltaFunc; ++index_DeltaFunc) {
                auto & temp_pos = Array_Pos_SmallDelta(index_DeltaFunc);
                if (temp_pos < _pos_Large) {
                    _array_info.emplace_back(Type_UpdateInfo_Small{index_DeltaFunc, temp_pos, _pos_Large});
                }

            }
        }

    };
}



