//
// Created by Angelina Chvirova on 07.03.2025.
//

#ifndef CUBESLY_NOISE_H
#define CUBESLY_NOISE_H

#include "TFile.h"
#include <iostream>
#include <TH1F.h>
#include <tuple>

class Noise {
public:
    explicit Noise(const std::string& inFileName);
    void SetTitlesHistNoise(const std::string& titleX, const std::string& titleY);
    std::pair<double, double> MaxPeak(unsigned int type, TH1F* hist);  /*Если гистограмма шума, то type = 0, если сигнала, то type = 1*/
    TH1F CutNoise(TH1F* hist_noise, TH1F* hist_LY);
    void WriteHistograms(unsigned int channel);
    ~Noise();

    std::map<int, TH1F>& GetFiberLYNoise() {
        return fiberLY_noise_;
    }
private:
    std::unique_ptr<TFile> file;
    const std::string inFileName_;
    std::map<int, TH1F> fiberLY_noise_;
};

#endif //CUBESLY_NOISE_H







