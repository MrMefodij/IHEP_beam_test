//
// Created by amefodev on 08.11.2024.
//

#ifndef CUBESLY_CALIBRATION_H
#define CUBESLY_CALIBRATION_H

#include <cmath>
#include <string>
#include <map>
// ch gain gain_error peak_1 peak_1_error mean_w/_background x-talk_w/_background mean_w/o_background x-talk_w/o_background
struct CalibResults{
    double _gain{NAN};
    double _gain_error{NAN};
    double _peak0{NAN};
    double _peak0_error{NAN};
    double _mean{NAN};
    double _xTalk{NAN};
};

class Calibration {
public:
    explicit Calibration(const std::string& inFileName, const int verbose=0);
    bool InitCalibByTXT_I();
    bool InitCalibByTXT_A();
    bool InitChannelFlagByTXT();
    const std::map<int, CalibResults>* GetCalibResultsI() const;
    const std::map<int, CalibResults>* GetCalibResultsA() const;
    const std::map<int, bool>* GetChannelFlag() const;

private:
    const std::string inFileName_;
    const int verbose_;
    std::string calibrationFileNameI_;
    std::string calibrationFileNameA_;
    bool calibration_I_Exist_;
    bool calibration_A_Exist_;
    bool channelFlag_Exist_;
    std::map<int, CalibResults> calibrationResultsI_;
    std::map<int, CalibResults> calibrationResultsA_;
    std::map<int, bool> channelFlag_;
};


#endif //CUBESLY_CALIBRATION_H
