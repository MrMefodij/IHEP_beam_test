//
// Created by Angelina Chvirova on 07.03.2025.
//

#include "Noise.h"
#include "FileHandler.h"
#include "defines.h"
#include "TF1.h"
#include <string>
#include <TFile.h>
#include <TSpectrum.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH1.h>
#include <filesystem>
#include <tuple>

Noise::Noise(const std::string& inFileName) : inFileName_(std::move(inFileName)) {
    file = std::make_unique<TFile>(inFileName.c_str(), "READ");

    if (file->IsZombie()) {
        std::cerr << "Не удалось открыть файл " << inFileName_ << std::endl;
        return;
    }

    TDirectory *dir = file.get();
    TIter next(dir->GetListOfKeys());
    TKey *key;

    while ((key = (TKey *) next())) {
        TObject *obj = key->ReadObj();

        if (!obj) continue;

        TDirectory *subdir = dynamic_cast<TDirectory *>(obj);
        if (subdir) {
            TIter subNext(subdir->GetListOfKeys());
            TKey *subKey;

            while ((subKey = (TKey *) subNext())) {
                TObject *subObj = subKey->ReadObj();
                if (!subObj || !subObj->InheritsFrom("TH1")) continue;

                std::string objName = subObj->GetName();
                TH1F *hist = dynamic_cast<TH1F *>(subObj);
                if (hist) {
                    TH1F newHist = *hist;
                    std::string name =subObj->GetName();
                    auto it = name.rfind(' ');
                    auto k = name.substr(it,name.size()-it);
                    int value = std::stoi(k);
                    if (objName.find("Noise") != std::string::npos) {
//                            std::cout << subObj->GetName() << " " <<AS_KV(newHist.GetEntries()) << " " << k <<"\n";
                    fiberLY_noise_[value] = newHist;
                    }
                }
            }
        }
    }
}

Noise::~Noise() {
    if (file) file->Close();
}

/*Если гистограмма шума, то type = 0, если сигнала, то type = 1*/
std::pair<double, double> Noise::MaxPeak(unsigned int type, TH1F* hist){
    std::unique_ptr<TSpectrum> fTSpectrum = std::make_unique<TSpectrum>(60);
    int nfound = 0;

    if (type == 0){
        nfound = fTSpectrum->Search(hist, 2, " ", 0.5);
    }
    if (type == 1){
        nfound = fTSpectrum->Search(hist, 1.0, " ", 0.1);
    }

    double first_peak = 0.0;
    double first_peak_Content = 0.0;
    if (nfound > 0) {
        double *xpeaks = fTSpectrum->GetPositionX();
        std::vector<double> peaks(xpeaks, xpeaks + nfound);
        std::sort(peaks.begin(), peaks.end());

        double first = peaks.front();
        int bin_first = hist->FindBin(first);
        double binCenter_first = hist->GetBinCenter(bin_first);
        double binContent_first  = hist->GetBinContent(bin_first);

        double second = peaks[1];
        int bin_second = hist->FindBin(second);
        double binContent_second = hist->GetBinContent(bin_second);

        if (binContent_first >= binContent_second){
            double binContent_left = hist->GetBinContent(bin_first - 1);
            double binContent_right = hist->GetBinContent(bin_first + 1);
            if (binContent_first > binContent_left && binContent_first > binContent_right){
                first_peak = binCenter_first;
                first_peak_Content = binContent_first;
            } else if (binContent_first <= binContent_left){
                first_peak = hist->GetBinCenter(bin_first - 1);
                first_peak_Content = binContent_left;
            }  else if (binContent_first <= binContent_right){
                first_peak = hist->GetBinCenter(bin_first + 1);
                first_peak_Content = binContent_right;
            }

        } else if (binContent_first < binContent_second){
            TF1 *fit = new TF1("fit", "gaus", first - 0.1, first + 0.1);
            hist->Fit("fit", "QR+");
            first_peak = fit->GetParameter(1);
            first_peak_Content = fit->GetParameter(0);
            }
        }
    return {first_peak, first_peak_Content};
}


TH1F Noise::CutNoise(TH1F* hist_noise, TH1F* hist_LY) {
    TH1F hist_LY_copy(*hist_LY);

    std::pair<double, double> f_noise = MaxPeak(0, hist_noise);
    double firstBinNoise = f_noise.first;
    double firstBinContentNoise = f_noise.second;

    std::pair<double, double> f_LY = MaxPeak(1, hist_LY);
    double firstBinLY = f_LY.first;
    double firstBinContentLY = f_LY.second;

//    std::cout << hist_noise->GetName() << " " << firstBinNoise << " " << firstBinContentNoise <<"\n";
//    std::cout << hist_LY_copy.GetName() << " " << firstBinLY <<  " " << firstBinContentLY <<"\n";

    double shift = firstBinLY - firstBinNoise;

//    if (shift != 0 ) {
//        shift = std::round(shift * 100) / 100;
//        std::cout << hist_LY_copy.GetName()  << " shift with noise:" << " " << shift <<"\n";
//    }

    int bins = hist_noise->GetNbinsX();
    TH1F hist_noise_shifted("hist_noise_shifted", hist_noise->GetTitle(), bins, 0, 200);
    for (int bin = 1; bin <= bins; ++bin) {
        double binContent = hist_noise->GetBinContent(bin);
        double binCenter = hist_noise->GetBinCenter(bin) + shift;
        int newBin = hist_noise_shifted.FindBin(binCenter);
        hist_noise_shifted.SetBinContent(newBin, binContent);
    }

    hist_LY_copy.Add(&hist_noise_shifted, - firstBinContentLY / firstBinContentNoise);
    hist_LY_copy.SetName(Form("%s_no_noise", hist_LY->GetName()));
    hist_LY_copy.SetTitle(Form("%s_no_noise", hist_LY->GetName()));

    for (int bin = 1; bin <= hist_LY_copy.GetNbinsX(); bin++) {
        double content = hist_LY_copy.GetBinContent(bin);
        if (content <= 0) {
            hist_LY_copy.SetBinContent(bin, 0);
        }
    }
    return hist_LY_copy;
}


void Noise::SetTitlesHistNoise(const std::string& titleX, const std::string& titleY) {
    std::array<std::map<int, TH1F>*, 1> mapHistograms = {&fiberLY_noise_};

    for (auto* mapPtr : mapHistograms) {
        if (mapPtr) {
            for (auto& [channel, hist] : *mapPtr) {
                hist.GetXaxis()->SetTitle(titleX.c_str());
                hist.GetYaxis()->SetTitle(titleY.c_str());
            }
        }
    }
}

void Noise::WriteHistograms(unsigned int channel) {
    std::vector<std::pair<TH1F&, std::string>> histograms = {
            {fiberLY_noise_[channel], "LY_noise"},
    };

    bool hasEntries = false;

    for (auto& [hist, suffix] : histograms) {
        if (hist.GetEntries() > 0) {
            std::string name = "fiber_" + std::to_string(channel) + "_" + suffix;
            hist.SetName(name.c_str());
            hist.SetTitle(name.c_str());
            hist.Write();
            hasEntries = true;
        }
    }
}

