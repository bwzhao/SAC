#include <iostream>
#include <vector>
#include "Config.h"
#include "Class_Calculate.h"
#include <sstream>
#include <iomanip>

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

    // bound for the omegas
    AC::Type_ValReal Min_Omega = std::stod(para_main[0 + index_StartPara]);
    AC::Type_ValReal Max_Omega = std::stod(para_main[1 + index_StartPara]);

    // Number of the grid
    int Num_DivideOmega = std::stoi(para_main[2 + index_StartPara]);

    // Number of the delta functions
    int Num_DeltaFunc = std::stoi(para_main[3 + index_StartPara]);

    // Number of grid in the measurement process
    int Num_GridBins = std::stoi(para_main[4 + index_StartPara]);

    // Where the files locates
    int index_StartFile = 6;
    auto str_dirInput = std::string("/Users/bowen/Dropbox/SandvikGroup/Projects/Current/Project_tJDynamics/DataAnalysis/2ndData/");
    auto str_dirOutput = std::string("/Users/bowen/Dropbox/SandvikGroup/Projects/Current/Project_tJDynamics/DataAnalysis/SpectralData/");
    double val_g = std::stod(para_main[index_StartFile + 0]);
    double val_h = std::stod(para_main[index_StartFile + 1]);
    int val_x = std::stoi(para_main[index_StartFile + 2]);
    int val_y = val_x;
    int val_q = std::stoi(para_main[index_StartFile + 3]);
    double val_beta = std::stod(para_main[index_StartFile + 4]);
    double val_delta=0.1;

    // Weight for the special leading deltq function
    int index_WeightLeading = 14;
    double WeightLeading = std::stod(para_main[0 + index_WeightLeading]);

    std::stringstream ss_fileName;
    ss_fileName << "_h" << std::fixed << std::setprecision(6) << val_h;
    ss_fileName << "_g" << std::fixed << std::setprecision(6) << val_g;
    ss_fileName << "_x" << val_x;
    ss_fileName << "_y" << val_y;
    ss_fileName << "_beta" << std::fixed << std::setprecision(6) << val_beta;
    ss_fileName << "_q" << val_q;
    ss_fileName << "_Delta" << std::fixed << std::setprecision(6) << val_delta;

    std::string file_val = str_dirInput + "val" + ss_fileName.str() + ".txt";
    std::string file_cov = str_dirInput + "cov" + ss_fileName.str() + ".txt";

    ss_fileName << "_Weight" << std::fixed << std::setprecision(6) << WeightLeading;
    std::string file_output_para = str_dirOutput + "para" + ss_fileName.str() + ".txt";
    std::string file_output = str_dirOutput + "data" + ss_fileName.str() + ".txt";

    // Parameters for measurement
    int index_StartMeaPara = 11;
    int Num_Sweeps = std::stoi(para_main[0 + index_StartMeaPara]);
    int Num_MeaBins = std::stoi(para_main[1 + index_StartMeaPara]);
    int Num_UpdateinSweep = std::stoi(para_main[2 + index_StartMeaPara]);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Start calculation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto LittleLion = AC::Class_Calculate(Min_Omega, Max_Omega, Num_DivideOmega,
                                          Num_DeltaFunc,
                                          file_val, file_cov, file_output_para, file_output,
                                          Num_GridBins,
                                          val_beta, WeightLeading);

    LittleLion.Anneal_Theta();
    LittleLion.CleanBin_Spectral();
    LittleLion.Setup_Measure();

    for (int index_Bins = 0; index_Bins != Num_MeaBins; ++index_Bins) {
        for (int index_Measure = 0; index_Measure != Num_Sweeps; ++index_Measure) {
            for (int index_UpdateinSweep = 0; index_UpdateinSweep != Num_UpdateinSweep; ++index_UpdateinSweep) {
                LittleLion.Update_SmallDeltaFunc(1, LittleLion.specFunc, LittleLion.array_G_tilde);
                LittleLion.Update_SmallDeltaFunc(2, LittleLion.specFunc, LittleLion.array_G_tilde);
                LittleLion.Update_LargeDeltaFunc(LittleLion.specFunc, LittleLion.array_G_tilde);
            }
            LittleLion.Measure_Spectral();
        }
        LittleLion.WriteBin_Spectral(index_Bins);
        LittleLion.Write_Para();
    }


    return 0;
}