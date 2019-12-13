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

    int index_StartFile = 5;
    std::string file_val = para_main[index_StartFile + 0];
    std::string file_cov = para_main[index_StartFile + 1];

    int index_StartModelPara = 7;
    AC::Type_ValReal Val_Beta = std::stod(para_main[0 + index_StartModelPara]);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Start calculation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto LittleLion = AC::Class_Calculate(Min_Omega, Max_Omega, Num_DivideOmega, Num_DeltaFunc, file_val, file_cov, 0);

    return 0;
}