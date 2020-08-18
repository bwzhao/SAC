#include <iostream>
#include <vector>
#include "Config.h"
#include "Class_Calculate.h"

int main(int argc, char *argv[]) {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Pre-calculation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //input parameters
    std::vector<std::string> para_main;
    for (decltype(argc) index_arg = 0; index_arg != argc; ++index_arg) {
        std::string temp_string(argv[index_arg]);
        para_main.push_back(temp_string);
    }
    int index_StartPara = 1;

    AC::Type_ValReal Min_Omega = std::stod(para_main[0 + index_StartPara]);
    AC::Type_ValReal Max_Omega = std::stod(para_main[1 + index_StartPara]);
    int Num_DivideOmega = std::stoi(para_main[2 + index_StartPara]);
    int Num_DeltaFunc = std::stoi(para_main[3 + index_StartPara]);
    int Num_Bins = std::stoi(para_main[4 + index_StartPara]);

    int index_StartFile = 6;
    std::string file_val = para_main[index_StartFile + 0];
    std::string file_cov = para_main[index_StartFile + 1];
    std::string file_output = para_main[index_StartFile + 2];

    int index_StartModelPara = 9;
    AC::Type_ValReal Val_Beta = std::stod(para_main[0 + index_StartModelPara]);

    int index_StartMeaPara = 10;
    int Num_Sweeps = std::stoi(para_main[0 + index_StartMeaPara]);
    int Num_MeaBins = std::stoi(para_main[1 + index_StartMeaPara]);
    int Num_UpdateinSweep = std::stoi(para_main[2 + index_StartMeaPara]);

    int index_WeightLeading = 13;
    double WeightLeading = std::stod(para_main[0 + index_WeightLeading]);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Start calculation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto LittleLion = AC::Class_Calculate(Min_Omega, Max_Omega, Num_DivideOmega, Num_DeltaFunc, file_val, file_cov, file_output,
                                          Num_Bins,
                                          Val_Beta);

    LittleLion.Anneal_Theta();
    LittleLion.CleanBin_Spectral();

    for (int index_Bins = 0; index_Bins != Num_MeaBins; ++index_Bins) {
        std::cout << index_Bins << std::endl;
        for (int index_Measure = 0; index_Measure != Num_Sweeps; ++index_Measure) {
            for (int index_UpdateinSweep = 0; index_UpdateinSweep != Num_UpdateinSweep; ++index_UpdateinSweep) {
                LittleLion.Update_DeltaFunc(1);
                LittleLion.Update_DeltaFunc(2);
            }
            LittleLion.Measure_Spectral();
        }
        LittleLion.WriteBin_Spectral();
    }

    return 0;
}