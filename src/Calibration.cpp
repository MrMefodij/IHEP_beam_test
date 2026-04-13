//
// Created by amefodev on 08.11.2024.
//

#include "Calibration.h"
#include "defines.h"
#include <fstream>
#include <iostream>

Calibration::Calibration(const std::string& inFileName, const int verbose) : inFileName_(inFileName), verbose_(verbose){
    size_t pos = inFileName.rfind('.');
    if (pos != std::string::npos) {
        calibrationFileNameI_ = inFileName_.substr(0, pos) + "_Calibration_I.txt";
        calibrationFileNameA_ = inFileName_.substr(0, pos) + "_Calibration_A.txt";
    }
    pos = calibrationFileNameI_.rfind('/');
    if (pos != std::string::npos) {
        calibrationFileNameI_ = "../calibration/I/" + calibrationFileNameI_.substr(pos, calibrationFileNameI_.size() - pos);
    }
    pos = calibrationFileNameA_.rfind('/');
    if (pos != std::string::npos) {
        calibrationFileNameA_ = "../calibration/A/" + calibrationFileNameA_.substr(pos, calibrationFileNameA_.size() - pos);
    }
    InitCalibByTXT_I();
    InitCalibByTXT_A();
    InitChannelFlagByTXT();
}

bool Calibration::InitCalibByTXT_I() {
    std::ifstream map(calibrationFileNameI_.c_str());
    if(!map.is_open()){
        std::cerr << calibrationFileNameI_ << " not found. Calibration is not available for that run.\n";
        calibration_I_Exist_ = false;
        return calibration_I_Exist_;
    }
    std::string header;
    std::getline(map, header);
    while (!map.eof()){
        int ch;
        map >> ch;
        if(std::floor(ch) == ch){
            map >> calibrationResultsI_[ch]._gain >> calibrationResultsI_[ch]._gain_error >> calibrationResultsI_[ch]._peak0 >>
            calibrationResultsI_[ch]._peak0_error >> calibrationResultsI_[ch]._mean >> calibrationResultsI_[ch]._xTalk;
            if (verbose_){
                std::cout << AS_KV(ch) << " " << calibrationResultsI_[ch]._gain << " " << calibrationResultsI_[ch]._gain_error << " " << calibrationResultsI_[ch]._peak0 << " " << calibrationResultsI_[ch]._peak0_error << std::endl;
            }
        }
    }
    calibration_I_Exist_ = true;
    return calibration_I_Exist_;
}

bool Calibration::InitCalibByTXT_A(){
    std::ifstream map(calibrationFileNameA_.c_str());
    if(!map.is_open()){
        std::cerr << calibrationFileNameA_ << " not found. Calibration is not available for that run.\n";
        calibration_A_Exist_ = false;
        return calibration_A_Exist_;
    }
    std::string header;
    std::getline(map, header);
    while (!map.eof()){
        int ch;
        map >> ch;
        if(std::floor(ch) == ch){
            calibrationResultsA_[ch]._mean =1; calibrationResultsA_[ch]._xTalk=1;
            map >> calibrationResultsA_[ch]._gain >> calibrationResultsA_[ch]._gain_error >> calibrationResultsA_[ch]._peak0 >>
            calibrationResultsA_[ch]._peak0_error >> calibrationResultsA_[ch]._mean >> calibrationResultsA_[ch]._xTalk;
            if (verbose_){
                std::cout << AS_KV(ch) << " " << calibrationResultsA_[ch]._gain << " " << calibrationResultsA_[ch]._gain_error << " " << calibrationResultsA_[ch]._peak0 << " " << calibrationResultsA_[ch]._peak0_error << std::endl;
            }
        }
    }
    calibration_A_Exist_ = true;
    return calibration_A_Exist_;
}

bool Calibration::InitChannelFlagByTXT() {
    std::ifstream flag(CHANNELS_FLAG);
    if(!flag.is_open()){
        std::cerr << CHANNELS_FLAG << " not found. Channels flag is not available.\n";
        channelFlag_Exist_ = false;
        return channelFlag_Exist_;
    }
    std::string header;
    std::getline(flag, header);
    while (!flag.eof()) {
        int ch;
        flag >> ch;
        if(std::floor(ch) == ch) {
            flag >> channelFlag_[ch];
            if (verbose_) {
                std::cout << AS_KV(ch) << " " << channelFlag_[ch] << std::endl;
            }
        }
    }
    channelFlag_Exist_ = true;
    return channelFlag_Exist_;
}

const std::map<int,CalibResults>* Calibration::GetCalibResultsI() const{
    if (calibration_I_Exist_)
        return &calibrationResultsI_;
    return nullptr;
}

const std::map<int,CalibResults>* Calibration::GetCalibResultsA() const{
    if (calibration_A_Exist_)
        return &calibrationResultsA_;
    return nullptr;
}

const std::map<int, bool>* Calibration::GetChannelFlag() const {
    if (channelFlag_Exist_)
        return &channelFlag_;
    return nullptr;
}