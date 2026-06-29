//
// Created by amefodev on 03.11.2024.
//

#include "HistogramManager.h"
#include <iostream>
#include <fstream>
#include <TF1.h>
#include <TSpectrum.h>
#include <memory>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <cmath>

HistogramManager::HistogramManager(int count){
    histograms_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "Cube_location_" + std::to_string(i);
        std::string title = "Cube_location_" + std::to_string(i);
        histograms_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }
}

TH2F* HistogramManager::getHistogram(size_t index) {
    if (index < histograms_.size()) {
        return histograms_[index].get();
    }
    return nullptr;
}

void HistogramManager::fillHistograms(unsigned int index, double x, double y){
    histograms_[index]->Fill(x,y);
}

void HistogramManager::SetBinContent(unsigned int index, double x, double y, double value){
    int binX = histograms_[index]->GetXaxis()->FindBin(x);
    int binY = histograms_[index]->GetYaxis()->FindBin(y);
    histograms_[index]->SetBinContent(binX,binY,value);
}


void HistogramManager::SetTitlesHist(const std::string& titleX, const std::string& titleY, const std::string& titleZ){
    for (auto& a : histograms_) {
        a->GetXaxis()->SetTitle(titleX.c_str());
        a->GetYaxis()->SetTitle(titleY.c_str());
        a->GetZaxis()->SetTitle(titleZ.c_str());
    }
}

void HistogramManager::WriteHistograms(TDirectory* dirName){
    dirName->cd();
    for (const auto& a : histograms_) {
        a->Write();
    }
}

void HistogramManager::FindHighDensityArea(int minEntries, std::ofstream& fileOut) {
    size_t cube = 0;
    fileOut << "cube x_low x_high y_low y_high\n";
    double globalMinX = std::numeric_limits<double>::max();
    double globalMaxX = std::numeric_limits<double>::lowest();
    double globalMinY = std::numeric_limits<double>::max();
    double globalMaxY = std::numeric_limits<double>::lowest();
    for (const auto& hist : histograms_) {

        if (!hist || hist->GetEntries() < minEntries) {
            std::cout << "[Skip] Cube " << cube << " has too few entries.\n";
            ++cube;
            continue;
        }

        double sum = 0.0;
        double currentThreshold = 0.0;
        int activeBins = 0, unactiveBins = 0;
        int nBinsX = hist->GetNbinsX();
        int nBinsY = hist->GetNbinsY();

        for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
                double content = hist->GetBinContent(ix, iy);
                if (content > 0) {
                    sum += content;
                    ++activeBins;
                } else {
                    ++unactiveBins;
                }
            }
        }

        if (activeBins == 0) { ++cube; continue; };
        double factor = unactiveBins / activeBins;

        double maxContent = hist->GetMaximum();
        double entries = hist->GetEntries();

        currentThreshold = maxContent * 0.15; // спец. кат по количеству событий
        // для ранов:
        //  88 - 92
        //  96 - 110        0.17
        // 125 - 128        0.17
        // 129 - 143        0.18
        // 150 - 168        0.15, 0.19 даёт много размеров около 8 - 8.5 см,
        // 188 - 193        0.19

        // Устанавливаем начальные значения для границ
        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();

        // Перебор всех бинов гистограммы
        for (int i = 1; i <= nBinsX; ++i) {
            for (int j = 1; j <= nBinsY; ++j) {
                double binContent = hist->GetBinContent(i, j);
                if (binContent > currentThreshold) {  // Проверка порога
                    double x = hist->GetXaxis()->GetBinCenter(i);
                    double y = hist->GetYaxis()->GetBinCenter(j);

                    // Обновление границ
                    if (x < minX) minX = x;
                    if (x > maxX) maxX = x;
                    if (y < minY) minY = y;
                    if (y > maxY) maxY = y;
                }
            }
        }

        // Вывод результатов
        if (minX != std::numeric_limits<double>::max()) {
            std::cout << "\nCube " << cube << ", dynamic threshold: " << currentThreshold << ", factor: " << factor << "\n";
//            std::cout << "  Ratio: " << ratio << ", mean: " << meanValue << "\n";

            std::cout << "The boundaries of a high-density area: \n" << AS_KV(cube) <<"\n";
            std::cout << "X: [" << minX << ", " << maxX << "]\n";
            std::cout << "Y: [" << minY << ", " << maxY << "]\n";
            std::cout << "Size of the area: "
                      << (maxX - minX) << " x " << (maxY - minY) << " mm"  << std::endl;
            fileOut << cube << " " << minX << " " << maxX  << " " << minY  << " " << maxY <<"\n";

            auto* box = new TBox(minX, minY, maxX, maxY);
            box->SetFillStyle(0);
            box->SetLineColor(kRed);
            box->SetLineWidth(2);
            hist->GetListOfFunctions()->Add(box);

            if (minX < globalMinX) globalMinX = minX;
            if (maxX > globalMaxX) globalMaxX = maxX;
            if (minY < globalMinY) globalMinY = minY;
            if (maxY > globalMaxY) globalMaxY = maxY;
        } else {
            std::cout << "There are no areas with a density higher than the specified threshold.\n";
        }
        ++cube;
    }
    fileOut << GLOBAL_CUBE << " " << globalMinX - 1  << " " << globalMaxX + 1 << " " << globalMinY - 1 << " " << globalMaxY + 1;
}

LyHistogramManager::LyHistogramManager(unsigned int count) : HistogramManager(count){
    histogramsLY_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "LY_" + std::to_string(i);
        std::string title = "LY_summ_" + std::to_string(i);
        histogramsLY_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    x_HistogramsLY_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_LY_" + std::to_string(i);
        std::string title = "x_LY_pe_" + std::to_string(i);
        x_HistogramsLY_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    y_HistogramsLY_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_LY_" + std::to_string(i);
        std::string title = "y_LY_pe_" + std::to_string(i);
        y_HistogramsLY_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    histogramsLY_normalized_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "LY_normalized" + std::to_string(i);
        std::string title = "LY_summ_normalized" + std::to_string(i);
        histogramsLY_normalized_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    x_HistogramsLY_normalized_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_LY_normalized" + std::to_string(i);
        std::string title = "x_LY_pe_normalized_" + std::to_string(i);
        x_HistogramsLY_normalized_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    y_HistogramsLY_normalized_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_LY_pe" + std::to_string(i);
        std::string title = "y_LY_pe_normalized_" + std::to_string(i);
        y_HistogramsLY_normalized_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    histograms_layer_.reserve(LAYERS_COUNT);
    for (size_t i = 0; i < LAYERS_COUNT; ++i) {
        std::string name = "Layer_" + std::to_string(i);
        std::string title = "Layer_" + std::to_string(i);
        histograms_layer_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    histogramsLY_normalized_layer_.reserve(LAYERS_COUNT);
    for (size_t i = 0; i < LAYERS_COUNT; ++i) {
        std::string name = "LY_normalized_layer_" + std::to_string(i);
        std::string title = "LY_summ_normalized_layer_" + std::to_string(i);
        histogramsLY_normalized_layer_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    x_HistogramsLY_normalized_layer_.reserve(LAYERS_COUNT);
    for (size_t i = 0; i < LAYERS_COUNT; ++i) {
        std::string name = "x_LY_normalized_layer_" + std::to_string(i);
        std::string title = "x_LY_pe_normalized_layer_" + std::to_string(i);
        x_HistogramsLY_normalized_layer_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    y_HistogramsLY_normalized_layer_.reserve(LAYERS_COUNT);
    for (size_t i = 0; i < LAYERS_COUNT; ++i) {
        std::string name = "y_LY_pe_layer_" + std::to_string(i);
        std::string title = "y_LY_pe_normalized_layer_" + std::to_string(i);
        y_HistogramsLY_normalized_layer_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    histograms_fiberLY_normalized_layer_.reserve(CHANNELS_NUMBER);
    for (size_t i = 0; i < CHANNELS_NUMBER; ++i) {
        const char* name = Form("fiber_%lu_LY_normalized_layer_pe", i);
        const char* title = Form("fiber_%lu_LY_normalized_layer_pe", i);
        histograms_fiberLY_normalized_layer_.emplace_back(std::make_unique<TH1F>(name, title, FIBER_BINNING));
    }

    x_fiberLY_normalized_position_.reserve(CHANNELS_NUMBER);
    for (size_t i = 0; i < CHANNELS_NUMBER; ++i) {
        const char* name = Form("x_fiber_%lu_LY_pe_Y_position", i);
        const char* title = Form("x_fiber_%lu_LY_pe_Y_position", i);
        x_fiberLY_normalized_position_.emplace_back(std::make_unique<TH1F>(name, title, CUBES_BINNING_Y));
    }

    y_fiberLY_normalized_position_.reserve(CHANNELS_NUMBER);
    for (size_t i = 0; i < CHANNELS_NUMBER; ++i) {
        const char* name = Form("y_fiber_%lu_LY_pe_X_position", i);
        const char* title = Form("y_fiber_%lu_LY_pe_X_position", i);
        y_fiberLY_normalized_position_.emplace_back(std::make_unique<TH1F>(name, title, CUBES_BINNING_X));
    }

    x_fiber_peak_comparison_.reserve(CHANNELS_NUMBER);
    for (size_t i = 0; i < CHANNELS_NUMBER; ++i) {
        const char* name = Form("x_fiber_%lu_peak_comparison_Y_position", i);
        const char* title = Form("x_fiber_%lu_peak_comparison_Y_position", i);
        x_fiber_peak_comparison_.emplace_back(std::make_unique<TH1F>(name, title, CUBES_BINNING_Y));
    }

    y_fiber_peak_comparison_.reserve(CHANNELS_NUMBER);
    for (size_t i = 0; i < CHANNELS_NUMBER; ++i) {
        const char* name = Form("y_fiber_%lu_peak_comparison_X_position", i);
        const char* title = Form("y_fiber_%lu_peak_comparison_X_position", i);
        y_fiber_peak_comparison_.emplace_back(std::make_unique<TH1F>(name, title, CUBES_BINNING_X));
    }

    x_Avr_Mean_fiberLY_normalized_position_.reserve(CHANNELS_NUMBER);
    for (size_t i = 0; i < CHANNELS_NUMBER; ++i) {
        const char* name = Form("x_fiber_%lu_Avr_Mean_LY_pe_Y_position", i);
        const char* title = Form("x_fiber_%lu_Avr_Mean_LY_pe_Y_position", i);
        x_Avr_Mean_fiberLY_normalized_position_.emplace_back(std::make_unique<TH1F>(name, title, CUBES_BINNING_Y));
    }

    y_Avr_Mean_fiberLY_normalized_position_.reserve(CHANNELS_NUMBER);
    for (size_t i = 0; i < CHANNELS_NUMBER; ++i) {
        const char* name = Form("y_fiber_%lu_Avr_Mean_LY_pe_X_position", i);
        const char* title = Form("y_fiber_%lu_Avr_Mean_LY_pe_X_position", i);
        y_Avr_Mean_fiberLY_normalized_position_.emplace_back(std::make_unique<TH1F>(name, title, CUBES_BINNING_X));
    }

    fiberLY_.reserve(CHANNELS_NUMBER);
    for (size_t i = 0; i < CHANNELS_NUMBER; ++i) {
        const char* name = Form("fiber_%lu_LY", i);
        const char* title = Form("fiber_%lu_LY", i);
        fiberLY_.emplace_back(std::make_unique<TH1F>(name, title, FIBER_BINNING));
    }

    histogramsLY_gaus_layer_.reserve(LAYERS_COUNT);
    for (size_t i = 0; i < LAYERS_COUNT; ++i) {
        std::string name = "Layer_Summ_LY_gaus" + std::to_string(i);
        std::string title = "Layer_Summ_LY_gaus" + std::to_string(i);
        histogramsLY_gaus_layer_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    x_histogramsLY_gaus_layer_.reserve(LAYERS_COUNT);
    for (size_t i = 0; i < LAYERS_COUNT; ++i) {
        std::string name = "x_Layer_Summ_LY_gaus" + std::to_string(i);
        std::string title = "x_Layer_Summ_LY_gaus" + std::to_string(i);
        x_histogramsLY_gaus_layer_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    y_histogramsLY_gaus_layer_.reserve(LAYERS_COUNT);
    for (size_t i = 0; i < LAYERS_COUNT; ++i) {
        std::string name = "y_Layer_Summ_LY_gaus" + std::to_string(i);
        std::string title = "y_Layer_Summ_LY_gaus" + std::to_string(i);
        y_histogramsLY_gaus_layer_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    histogramsAverageLY_gaus_layer_ = std::make_unique<TH2F>("histogramsAverageLY_gaus_layer", "Gaus Mean LY X+Y", CUBES_BINNING);
    x_histogramsAverageLY_gaus_layer_ = std::make_unique<TH2F>("x_histogramsAverageLY_gaus_layer", "Gaus Mean LY X", CUBES_BINNING);
    y_histogramsAverageLY_gaus_layer_ = std::make_unique<TH2F>("y_histogramsAverageLY_gaus_layer", "Gaus Mean LY Y", CUBES_BINNING);

    histogramsLY_normalized_layer_slice_X_.reserve(CUBES_SLICE_X);
    histogramsLY_normalized_layer_slice_Y_.reserve(CUBES_SLICE_Y);
    y_HistogramsLY_normalized_layer_slice_.reserve(CUBES_SLICE_Y);
    x_HistogramsLY_normalized_layer_slice_.reserve(CUBES_SLICE_X);

    x_CubeLY_I.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_LY_I" + std::to_string(i);
        std::string title = "x_LY_I_pe_" + std::to_string(i);
        x_CubeLY_I.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), FIBER_BINNING));
    }

    y_CubeLY_I.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_LY_I" + std::to_string(i);
        std::string title = "y_LY_I_pe_" + std::to_string(i);
        y_CubeLY_I.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), FIBER_BINNING));
    }

    x_CubeLY_A.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_LY_A" + std::to_string(i);
        std::string title = "x_LY_A_pe_" + std::to_string(i);
        x_CubeLY_A.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), FIBER_BINNING));
    }

    y_CubeLY_A.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_LY_A" + std::to_string(i);
        std::string title = "y_LY_A_pe_" + std::to_string(i);
        y_CubeLY_A.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), FIBER_BINNING));
    }

    x_CubeLY_A_recon.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_LY_A_reconstructed" + std::to_string(i);
        std::string title = "x_LY_A_pe_reconstructed_" + std::to_string(i);
        x_CubeLY_A_recon.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), FIBER_BINNING));
    }

    y_CubeLY_A_recon.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_LY_A_reconstructed" + std::to_string(i);
        std::string title = "y_LY_A_pe_reconstructed_" + std::to_string(i);
        y_CubeLY_A_recon.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), FIBER_BINNING));
    }

    averageLY_ = std::make_unique<TH1F>("Average_LY", "Average_LY", FIBER_BINNING);
    averageLY_X_ = std::make_unique<TH1F>("Average_X_LY", "Average_X_LY", FIBER_BINNING);
    averageLY_Y_ = std::make_unique<TH1F>("Average_Y_LY", "Average_Y_LY", FIBER_BINNING);
    averageTime_ = std::make_unique<TH1F>("Average_Time", "Average_Time", 600, 0 ,600);
    averageCubeResponse_ =  std::make_unique<TH2F>("Average_CubeResponse", "Average_CubeResponse", CUBES_BINNING);

//    pixelsAverageLY_X_Y_ = std::make_unique<TH1F>("Pixel_Average_LY_X_Y", "Pixel_Average_LY_X_Y", CUBES_BINNING_Y);
    pixelsAverageLY_X_ = std::make_unique<TH1F>("Pixel_Average_LY_X", "Pixel_Average_LY_X", CUBES_BINNING_X); // 1D
    pixelsAverageLY_Y_ = std::make_unique<TH1F>("Pixel_Average_LY_Y", "Pixel_Average_LY_Y", CUBES_BINNING_Y);

    averageLY_->GetXaxis()->SetTitle("LY(X or Y)[p.e.]");
    averageLY_->GetYaxis()->SetTitle("Number");

    averageLY_X_->GetXaxis()->SetTitle("LY(X)[p.e.]");
    averageLY_X_->GetYaxis()->SetTitle("Number");

    averageLY_Y_->GetXaxis()->SetTitle("LY(Y)[p.e.]");
    averageLY_Y_->GetYaxis()->SetTitle("Number");

    averageTime_->GetXaxis()->SetTitle("t[8ns]");
    averageTime_->GetYaxis()->SetTitle("Number");
}

void LyHistogramManager::SetBinContentLY(unsigned int index, double x, double y, double x_value, double y_value){
    int binXly = histogramsLY_[index]->GetXaxis()->FindBin(x);
    int binYly = histogramsLY_[index]->GetYaxis()->FindBin(y);
    auto currentValue = histogramsLY_[index]->GetBinContent(binXly,binYly);

    int x_binXly = x_HistogramsLY_[index]->GetXaxis()->FindBin(x);
    int x_binYly = x_HistogramsLY_[index]->GetYaxis()->FindBin(y);
    auto x_currentValue = x_HistogramsLY_[index]->GetBinContent(x_binXly,x_binYly);

    int y_binXly = y_HistogramsLY_[index]->GetXaxis()->FindBin(x);
    int y_binYly = y_HistogramsLY_[index]->GetYaxis()->FindBin(y);
    auto y_currentValue = y_HistogramsLY_[index]->GetBinContent(y_binXly,y_binYly);

    int binX = histograms_[index]->GetXaxis()->FindBin(x);
    int binY = histograms_[index]->GetYaxis()->FindBin(y);
    auto currentCount = histograms_[index]->GetBinContent(binX,binY);

    histograms_[index]->Fill(x,y);
    double newAvr = (currentValue * currentCount + x_value + y_value)/(currentCount + 1);
    histogramsLY_[index]->SetBinContent(binXly,binYly, newAvr);

    double x_newAvr = (x_currentValue * currentCount + x_value)/(currentCount + 1);
    x_HistogramsLY_[index]->SetBinContent(binXly,binYly, x_newAvr);

    double y_newAvr = (y_currentValue * currentCount + y_value)/(currentCount + 1);
    y_HistogramsLY_[index]->SetBinContent(binXly,binYly, y_newAvr);
    averageLY_->Fill(x_value);
    averageLY_X_->Fill(x_value);
    averageLY_->Fill(y_value);
    averageLY_Y_->Fill(y_value);
}

void LyHistogramManager::SetBinContentNormalizedLY(unsigned int index, double x, double y, double x_value, double y_value, double z_value){
    int binXly = histogramsLY_normalized_[index]->GetXaxis()->FindBin(x);
    int binYly = histogramsLY_normalized_[index]->GetYaxis()->FindBin(y);
    auto currentValue_norm = histogramsLY_normalized_[index]->GetBinContent(binXly,binYly);

    int x_binXly = x_HistogramsLY_normalized_[index]->GetXaxis()->FindBin(x);
    int x_binYly = x_HistogramsLY_normalized_[index]->GetYaxis()->FindBin(y);
    auto x_currentValue_norm = x_HistogramsLY_normalized_[index]->GetBinContent(x_binXly,x_binYly);

    int y_binXly = y_HistogramsLY_normalized_[index]->GetXaxis()->FindBin(x);
    int y_binYly = y_HistogramsLY_normalized_[index]->GetYaxis()->FindBin(y);
    auto y_currentValue_norm = y_HistogramsLY_normalized_[index]->GetBinContent(y_binXly,y_binYly);

    int binX = histograms_[index]->GetXaxis()->FindBin(x);
    int binY = histograms_[index]->GetYaxis()->FindBin(y);
    auto currentCount = histograms_[index]->GetBinContent(binX,binY);

    histograms_[index]->Fill(x,y);
    double newAvr = (currentValue_norm * currentCount + x_value + y_value + z_value)/(currentCount + 1);
    histogramsLY_normalized_[index]->SetBinContent(binXly,binYly, newAvr);

    double x_newAvr_norm = (x_currentValue_norm * currentCount + x_value)/(currentCount + 1);
    x_HistogramsLY_normalized_[index]->SetBinContent(binXly,binYly, x_newAvr_norm);

    double y_newAvr_norm = (y_currentValue_norm * currentCount + y_value)/(currentCount + 1);
    y_HistogramsLY_normalized_[index]->SetBinContent(binXly,binYly, y_newAvr_norm);
}

void LyHistogramManager::SetBinContentNormalizedLY(unsigned int index, double x, double y, double x_value, double y_value){
    int binXly = histogramsLY_normalized_[index]->GetXaxis()->FindBin(x);
    int binYly = histogramsLY_normalized_[index]->GetYaxis()->FindBin(y);
    auto currentValue_norm = histogramsLY_normalized_[index]->GetBinContent(binXly,binYly);

    auto x_currentValue_norm = x_HistogramsLY_normalized_[index]->GetBinContent(binXly,binYly);
    auto y_currentValue_norm = y_HistogramsLY_normalized_[index]->GetBinContent(binXly,binYly);

    int binX = histograms_[index]->GetXaxis()->FindBin(x);
    int binY = histograms_[index]->GetYaxis()->FindBin(y);
    auto currentCount = histograms_[index]->GetBinContent(binX,binY);

    histograms_[index]->Fill(x,y);
    double newAvr = (currentValue_norm * currentCount + x_value + y_value)/(currentCount + 1);
    histogramsLY_normalized_[index]->SetBinContent(binXly,binYly, newAvr);

    double x_newAvr_norm = (x_currentValue_norm * currentCount + x_value)/(currentCount + 1);
    x_HistogramsLY_normalized_[index]->SetBinContent(binXly,binYly, x_newAvr_norm);

    double y_newAvr_norm = (y_currentValue_norm * currentCount + y_value)/(currentCount + 1);
    y_HistogramsLY_normalized_[index]->SetBinContent(binXly,binYly, y_newAvr_norm);
}

void LyHistogramManager::SetBinContentNormalizedLY_layer(unsigned int index, double x, double y, double x_value, double y_value) {
    int binXly = histogramsLY_normalized_layer_[index]->GetXaxis()->FindBin(x);
    int binYly = histogramsLY_normalized_layer_[index]->GetYaxis()->FindBin(y);
    auto currentValue_norm = histogramsLY_normalized_layer_[index]->GetBinContent(binXly,binYly);

    auto x_currentValue_norm = x_HistogramsLY_normalized_layer_[index]->GetBinContent(binXly,binYly);
    auto y_currentValue_norm = y_HistogramsLY_normalized_layer_[index]->GetBinContent(binXly,binYly);
    auto currentCount = histograms_layer_[index]->GetBinContent(binXly,binYly);

    histograms_layer_[index]->Fill(x,y);
    double newAvr = (currentValue_norm * currentCount + x_value + y_value)/(currentCount + 1);
    histogramsLY_normalized_layer_[index]->SetBinContent(binXly,binYly, newAvr);

    double x_newAvr_norm = (x_currentValue_norm * currentCount + x_value)/(currentCount + 1);
    x_HistogramsLY_normalized_layer_[index]->SetBinContent(binXly,binYly, x_newAvr_norm);

    double y_newAvr_norm = (y_currentValue_norm * currentCount + y_value)/(currentCount + 1);
    y_HistogramsLY_normalized_layer_[index]->SetBinContent(binXly,binYly, y_newAvr_norm);

    pixelsDataSummLY_[index][binXly][binYly].push_back(x_value + y_value);
    pixelsDataXLY_[index][binXly][binYly].push_back(x_value);
    pixelsDataYLY_[index][binXly][binYly].push_back(y_value);
}

void LyHistogramManager::SetBinContentNormalizedLY_layer(unsigned int index, double x, double y, double x_value, double y_value, double z_value) {
    int binXly = histogramsLY_normalized_layer_[index]->GetXaxis()->FindBin(x);
    int binYly = histogramsLY_normalized_layer_[index]->GetYaxis()->FindBin(y);
    auto currentValue_norm = histogramsLY_normalized_layer_[index]->GetBinContent(binXly,binYly);

    auto x_currentValue_norm = x_HistogramsLY_normalized_layer_[index]->GetBinContent(binXly,binYly);
    auto y_currentValue_norm = y_HistogramsLY_normalized_layer_[index]->GetBinContent(binXly,binYly);

    int binX = histograms_layer_[index]->GetXaxis()->FindBin(x);
    int binY = histograms_layer_[index]->GetYaxis()->FindBin(y);
    auto currentCount = histograms_layer_[index]->GetBinContent(binX,binY);

    histograms_layer_[index]->Fill(x,y);
    double newAvr = (currentValue_norm * currentCount + x_value + y_value + z_value)/(currentCount + 1);
    histogramsLY_normalized_layer_[index]->SetBinContent(binXly,binYly, newAvr);

    double x_newAvr_norm = (x_currentValue_norm * currentCount + x_value)/(currentCount + 1);
    x_HistogramsLY_normalized_layer_[index]->SetBinContent(binXly,binYly, x_newAvr_norm);

    double y_newAvr_norm = (y_currentValue_norm * currentCount + y_value)/(currentCount + 1);
    y_HistogramsLY_normalized_layer_[index]->SetBinContent(binXly,binYly, y_newAvr_norm);

    pixelsDataSummLY_[index][binXly][binYly].push_back(x_value + y_value + z_value);
    pixelsDataXLY_[index][binXly][binYly].push_back(x_value);
    pixelsDataYLY_[index][binXly][binYly].push_back(y_value);
}


void LyHistogramManager::SetFiberLY_layer(double x_fiber, double y_fiber, double x_position, double y_position,  double x_ly, double y_ly){
    histograms_fiberLY_normalized_layer_[static_cast<size_t>(x_fiber)]->Fill(x_ly);
    histograms_fiberLY_normalized_layer_[static_cast<size_t>(y_fiber)]->Fill(y_ly);

    int binXly = histogramsLY_normalized_layer_[0]->GetXaxis()->FindBin(x_position);
    int binYly = histogramsLY_normalized_layer_[0]->GetYaxis()->FindBin(y_position);

    x_fiber_pixel_LY_[static_cast<size_t>(x_fiber)][binYly].push_back(x_ly);
    y_fiber_pixel_LY_[static_cast<size_t>(y_fiber)][binXly].push_back(y_ly);
}

void LyHistogramManager::SetTitlesHistLY(const std::string& titleX, const std::string& titleY, const std::string& titleZ){
    for (auto& a : histogramsLY_) {
        a->GetXaxis()->SetTitle(titleX.c_str());
        a->GetYaxis()->SetTitle(titleY.c_str());
        a->GetZaxis()->SetTitle(titleZ.c_str());
    }
    for (auto& a : x_HistogramsLY_) {
        a->GetXaxis()->SetTitle(titleX.c_str());
        a->GetYaxis()->SetTitle(titleY.c_str());
        a->GetZaxis()->SetTitle(("X_" + titleZ).c_str());
    }
    for (auto& a : y_HistogramsLY_) {
        a->GetXaxis()->SetTitle(titleX.c_str());
        a->GetYaxis()->SetTitle(titleY.c_str());
        a->GetZaxis()->SetTitle(("Y_" + titleZ).c_str());
    }
    for (auto& a : histogramsLY_normalized_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle((titleZ).c_str());
    }
    for (auto& a : x_HistogramsLY_normalized_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle(("X_" + titleZ).c_str());
    }
    for (auto& a : y_HistogramsLY_normalized_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle(("Y_" + titleZ).c_str());
    }
    for (auto& a : histogramsLY_normalized_layer_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle((titleZ).c_str());
    }
    for (auto& a : x_HistogramsLY_normalized_layer_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle(("X_" + titleZ).c_str());
    }
    for (auto& a : y_HistogramsLY_normalized_layer_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle(("Y_" + titleZ).c_str());
    }

    for (auto& a : histograms_fiberLY_normalized_layer_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
    }

    for (auto& a : x_fiberLY_normalized_position_) {
        a->GetXaxis()->SetTitle("Y[mm]");
        a->GetYaxis()->SetTitle("LY[p.e.]");
    }

    for (auto& a : x_fiber_peak_comparison_) {
        a->GetXaxis()->SetTitle("Y[mm]");
        a->GetYaxis()->SetTitle("peak_1 (LY < 18 pe) / peak_2 (LY >= 18 pe)");
    }

    for (auto& a : y_fiber_peak_comparison_) {
        a->GetXaxis()->SetTitle("X[mm]");
        a->GetYaxis()->SetTitle("peak_1 (LY < 18 pe) / peak_2 (LY >= 18 pe)");
    }

    for (auto& a : fiberLY_) {
        a->GetXaxis()->SetTitle("LY[p.e.]");
        a->GetYaxis()->SetTitle("Number");
    }

    for (auto& a : y_fiberLY_normalized_position_) {
        a->GetXaxis()->SetTitle("X[mm]");
        a->GetYaxis()->SetTitle("LY[p.e.]");
    }

    for (auto& a : x_Avr_Mean_fiberLY_normalized_position_) {
        a->GetXaxis()->SetTitle("Y[mm]");
        a->GetYaxis()->SetTitle("Avr LY[p.e.]");
    }

    for (auto& a : y_Avr_Mean_fiberLY_normalized_position_) {
        a->GetXaxis()->SetTitle("X[mm]");
        a->GetYaxis()->SetTitle("Avr LY[p.e.]");
    }

    for (auto& a : histogramsLY_gaus_layer_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle(("X+Y fibres" + titleZ).c_str());
    }

    for (auto& a : x_histogramsLY_gaus_layer_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle(("X fibre " + titleZ).c_str());
    }

    for (auto& a : y_histogramsLY_gaus_layer_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle(("Y fibre " + titleZ).c_str());
    }

    histogramsAverageLY_gaus_layer_->GetXaxis()->SetTitle((titleX).c_str());
    histogramsAverageLY_gaus_layer_->GetYaxis()->SetTitle((titleY).c_str());
    histogramsAverageLY_gaus_layer_->GetZaxis()->SetTitle(("X+Y_" + titleZ).c_str());

    x_histogramsAverageLY_gaus_layer_->GetXaxis()->SetTitle((titleX).c_str());
    x_histogramsAverageLY_gaus_layer_->GetYaxis()->SetTitle((titleY).c_str());
    x_histogramsAverageLY_gaus_layer_->GetZaxis()->SetTitle(("X_" + titleZ).c_str());

    y_histogramsAverageLY_gaus_layer_->GetXaxis()->SetTitle((titleX).c_str());
    y_histogramsAverageLY_gaus_layer_->GetYaxis()->SetTitle((titleY).c_str());
    y_histogramsAverageLY_gaus_layer_->GetZaxis()->SetTitle(("Y_" + titleZ).c_str());
}

void LyHistogramManager::SetTitlesHistLY_Fiber_I(const std::string& titleX, const std::string& titleY){
    for (int i = 0; i < CUBES_COUNT; ++i) {
        x_CubeLY_I[i]->GetYaxis()->SetTitle(titleY.c_str());
        y_CubeLY_I[i]->GetYaxis()->SetTitle(titleY.c_str());
        x_CubeLY_I[i]->GetXaxis()->SetTitle(titleX.c_str());
        y_CubeLY_I[i]->GetXaxis()->SetTitle(titleX.c_str());
    }
}
void LyHistogramManager::SetTitlesHistLY_Fiber_A(const std::string& titleX, const std::string& titleY){
    for (int i = 0; i < CUBES_COUNT; ++i) {
        x_CubeLY_A[i]->GetYaxis()->SetTitle(titleY.c_str());
        y_CubeLY_A[i]->GetYaxis()->SetTitle(titleY.c_str());
        x_CubeLY_A[i]->GetXaxis()->SetTitle(titleX.c_str());
        y_CubeLY_A[i]->GetXaxis()->SetTitle(titleX.c_str());
    }
}

void LyHistogramManager::SetTitlesHistLY_Fiber_A_reconstructed(const std::string& titleX, const std::string& titleY){
    for (int i = 0; i < CUBES_COUNT; ++i) {
        x_CubeLY_A_recon[i]->GetYaxis()->SetTitle(titleY.c_str());
        y_CubeLY_A_recon[i]->GetYaxis()->SetTitle(titleY.c_str());
        x_CubeLY_A_recon[i]->GetXaxis()->SetTitle(titleX.c_str());
        y_CubeLY_A_recon[i]->GetXaxis()->SetTitle(titleX.c_str());
    }
}

void LyHistogramManager::SetFiberLY_I(unsigned int index, double x, double y){
    x_CubeLY_I[index]->Fill(x);
    y_CubeLY_I[index]->Fill(y);
}

void LyHistogramManager::SetFiberLY_A(unsigned int index, double x, double y){
    x_CubeLY_A[index]->Fill(x);
    y_CubeLY_A[index]->Fill(y);
}

void LyHistogramManager::SetFiberLY_A_recon(unsigned int index, double x, double y){
    x_CubeLY_A_recon[index]->Fill(x);
    y_CubeLY_A_recon[index]->Fill(y);
}

void LyHistogramManager::FillFiberLY(const unsigned int ch_x, const unsigned int ch_y, const double ly_x, const double ly_y) {
//void LyHistogramManager::FillFiberLY(const unsigned int ch_x, const unsigned int ch_y, const unsigned int ch_z, const double ly_x, const double ly_y, const double ly_z) { // для z волокон
    fiberLY_[ch_x]->Fill(ly_x);
    fiberLY_[ch_y]->Fill(ly_y);
//    fiberLY_[ch_z]->Fill(ly_z);
}

void LyHistogramManager::FillEventTime(double value){
    averageTime_->Fill(value);
}

void LyHistogramManager::SetBinContentAvrCubeResponse(double x, double y, double value){

}

void LyHistogramManager::FillAveragePlots(const double x, const double y, const double x_ly, const double y_ly) { // может тут проблема
    const int binX = histogramsAverageLY_gaus_layer_->GetXaxis()->FindBin(x);
    const int binY = histogramsAverageLY_gaus_layer_->GetYaxis()->FindBin(y);

    pixelsAverageSummLY_[binX][binY].push_back(x_ly+y_ly);
    pixelsAverageXLY_[binX][binY].push_back(x_ly);
    pixelsAverageYLY_[binX][binY].push_back(y_ly);
}

void LyHistogramManager::NormalizeCubes() const { // нормировка по кубам
    auto fTSpectrum = new TSpectrum(15);
    int nfound = 0;
    TF1* fit;

    std::ofstream fileOut(CUBES_NORMALIZATION);
    fileOut << "cube side norm ( summ: ";
    double ly_norm = 0;
    double ly_norm_y = 0;
    double ly_norm_x = 0;

    nfound = fTSpectrum->Search(averageLY_.get(),10,"",0.5);
    double *LyPeaks = fTSpectrum->GetPositionX();
    for (uint peak = 0; peak < nfound; ++peak) {
        if (LyPeaks[peak] > 15 && LyPeaks[peak] < 75) {
            fit = new TF1("gaus","gaus", LyPeaks[peak] - 7, LyPeaks[peak] + 7);
            averageLY_->Fit(fit, "qr+");
            ly_norm = fit->GetParameter(1);
            fileOut << ly_norm <<", x: ";
        }
    }

    nfound = fTSpectrum->Search(averageLY_X_.get(),10,"",0.5);
    LyPeaks = fTSpectrum->GetPositionX();
    for (uint peak = 0; peak < nfound; ++peak) {
        if (LyPeaks[peak] > 15 && LyPeaks[peak] < 75) {
            fit = new TF1("gaus","gaus", LyPeaks[peak] - 7, LyPeaks[peak] + 7);
            averageLY_X_->Fit(fit, "qr+");
            ly_norm_x = fit->GetParameter(1);
            fileOut << ly_norm_x <<"; y: ";
        }
    }

    nfound = fTSpectrum->Search(averageLY_Y_.get(),10,"",0.5);
    LyPeaks = fTSpectrum->GetPositionX();
    for (uint peak = 0; peak < nfound; ++peak) {
        if (LyPeaks[peak] > 15 && LyPeaks[peak] < 75) {
            fit = new TF1("gaus","gaus", LyPeaks[peak] - 7, LyPeaks[peak] + 7);
            averageLY_Y_->Fit(fit, "qr+");
            ly_norm_y = fit->GetParameter(1);
            fileOut << ly_norm_y <<")\n";
        }
    }

    for (int i = 0; i < CUBES_COUNT; ++i) {
        nfound = fTSpectrum->Search(x_CubeLY_A_recon[i].get(),10,"",0.05);
        if(nfound > 0) {
            double *xpeaks = fTSpectrum->GetPositionX();
            for (uint peak = 0; peak < nfound; ++peak) {
                if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                    fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                    x_CubeLY_A_recon[i]->Fit(fit, "qr+");
                    fileOut << i << " x " << fit->GetParameter(1)/ly_norm << std::endl;
                    break;
                }
            }
        }

        nfound = fTSpectrum->Search(y_CubeLY_A_recon[i].get(),10,"",0.05);
        if(nfound > 0) {
            double *xpeaks = fTSpectrum->GetPositionX();
            for (uint peak = 0; peak < nfound; ++peak) {
                if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                    fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                    y_CubeLY_A_recon[i]->Fit(fit, "qr+");
                    fileOut << i << " y " << fit->GetParameter(1)/ly_norm << std::endl;
                    break;
                }
            }
        }
    }
    fileOut.close();
}

void LyHistogramManager::NormalizeFibers() const{ // нормировка по волокнам
    auto fTSpectrum = new TSpectrum(15);
    int nfound = 0;
    TF1* fit;

    std::ofstream fileOut(FIBERS_NORMALIZATION);
    fileOut << "fiber norm (avr: ";
    double ly_norm = 0;

    nfound = fTSpectrum->Search(averageLY_.get(),10,"",0.5);
    double *LyPeaks = fTSpectrum->GetPositionX();
    for (uint peak = 0; peak < nfound; ++peak) {
        if (LyPeaks[peak] > 15 && LyPeaks[peak] < 75) {
            fit = new TF1("gaus","gaus", LyPeaks[peak] - 7 , LyPeaks[peak] + 7);
            averageLY_->Fit(fit, "qr+");
            ly_norm = fit->GetParameter(1);
            fileOut << ly_norm <<")\n";
        }
    }

    for (int ch = 0; ch < CHANNELS_NUMBER; ++ch) {
        if (fiberLY_[ch]->GetEntries() > 1000) {
            nfound = fTSpectrum->Search(fiberLY_[ch].get(),10,"",0.05);
            if(nfound > 0) {
                double *xpeaks = fTSpectrum->GetPositionX();
                for (uint peak = 0; peak < nfound; ++peak) {
                    if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                        fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                        fiberLY_[ch]->Fit(fit, "qr+");
                        fileOut << ch << " " << fit->GetParameter(1)/ly_norm << std::endl;
                        break;
                    }
                }
            }
        } else {
            fileOut << ch << " " << 1 << std::endl;
        }
    }
}

void LyHistogramManager::WriteHistograms(TDirectory* dirName){
    auto count = dirName->mkdir("count");
    count->cd();
    for (const auto& a : histograms_) {
        a->Write();
    }

    auto avrLY = dirName->mkdir("avrLY");
    avrLY->cd();
    for (const auto& a : histogramsLY_) {
        a->Write();
    }

    auto avrLY_x = dirName->mkdir("x_avrLY");
    avrLY_x->cd();
    for (const auto& a : x_HistogramsLY_) {
        a->Write();
    }

    auto avrLY_y = dirName->mkdir("y_avrLY");
    avrLY_y->cd();
    for (const auto& a : y_HistogramsLY_) {
        a->Write();
    }

    auto avrLY_norm = dirName->mkdir("avrLY_norm");
    avrLY_norm->cd();
    for (const auto& a : histogramsLY_normalized_) {
        a->Write();
    }

    auto avrLY_x_norm = dirName->mkdir("x_avrLY_norm");
    avrLY_x_norm->cd();
    for (const auto& a : x_HistogramsLY_normalized_) {
        a->Write();
    }

    auto avrLY_y_norm = dirName->mkdir("y_avrLY_norm");
    avrLY_y_norm->cd();
    for (const auto& a : y_HistogramsLY_normalized_) {
        a->Write();
    }

    auto avrLY_layer = dirName->mkdir("avrLY_layer");
    avrLY_layer->cd();
    for (const auto& a : histograms_layer_) {
        a->GetXaxis()->SetTitle("X[mm]");
        a->GetYaxis()->SetTitle("Y[mm]");
        a->GetZaxis()->SetTitle("Number");
        a->Write();
    }
    // for (const auto& a : histogramsLY_normalized_layer_) {
    for (int j = 0; j < 3; ++j){
        avrLY_layer->cd();
        histogramsLY_normalized_layer_.at(j)->Write();
        auto avrLY_layer_slices_X = avrLY_layer->mkdir(Form("slices_X_LY_sum_layer_%d", j));
        avrLY_layer_slices_X->cd();
        for (int x = 0; x < CUBES_SLICE_X - 0; ++x) {
            histogramsLY_normalized_layer_slice_X_[x] = std::unique_ptr<TH1D>(histogramsLY_normalized_layer_.at(j)->ProjectionY(Form("slice_X_%d_layer_%d", x, j), x, x));
            if (histogramsLY_normalized_layer_slice_X_[x]->GetEntries() > 10) {
                histogramsLY_normalized_layer_slice_X_[x]->Write();
            }
        }
        auto avrLY_layer_slices_Y = avrLY_layer->mkdir(Form("slices_Y_LY_sum_layer_%d", j));
        avrLY_layer_slices_Y->cd();
        for (int y = 0; y < CUBES_SLICE_Y - 0; ++y) {
            histogramsLY_normalized_layer_slice_Y_[y] = std::unique_ptr<TH1D>(histogramsLY_normalized_layer_.at(j)->ProjectionX(Form("slice_Y_%d_layer_%d", y, j), y, y));
            if (histogramsLY_normalized_layer_slice_Y_[y]->GetEntries() > 10) {
                histogramsLY_normalized_layer_slice_Y_[y]->Write();
            }
        }
    }
    // for (const auto& a : x_HistogramsLY_normalized_layer_) {
    for (int j = 0; j < 3; ++j){
        avrLY_layer->cd();
        x_HistogramsLY_normalized_layer_.at(j)->Write();
        auto avrLY_layer_slices_X = avrLY_layer->mkdir(Form("avr_slices_X_LY_sum_layer_%d", j));
        avrLY_layer_slices_X->cd();
        for (int y = 0; y < CUBES_SLICE_Y - 0; ++y) {
            x_HistogramsLY_normalized_layer_slice_[y] = std::unique_ptr<TH1D>(x_HistogramsLY_normalized_layer_.at(j)->ProjectionX(Form("x_slice_Y_%d_layer_%d", y, j), y, y));
            if (x_HistogramsLY_normalized_layer_slice_[y]->GetEntries() > 10) {
                x_HistogramsLY_normalized_layer_slice_[y]->Write();
            }
        }
    }
    // for (const auto& a : y_HistogramsLY_normalized_layer_) {
    for (int j = 0; j < 3; ++j){
        avrLY_layer->cd();
        y_HistogramsLY_normalized_layer_.at(j)->Write();
        auto avrLY_layer_slices_Y = avrLY_layer->mkdir(Form("avr_slices_Y_LY_sum_layer_%d", j));
        avrLY_layer_slices_Y->cd();
        for (int x = 0; x < CUBES_SLICE_Y - 0; ++x) {
            y_HistogramsLY_normalized_layer_slice_[x] = std::unique_ptr<TH1D>(y_HistogramsLY_normalized_layer_.at(j)->ProjectionX(Form("y_slice_X_%d_layer_%d", x, j), x, x));
            if (y_HistogramsLY_normalized_layer_slice_[x]->GetEntries() > 10) {
                y_HistogramsLY_normalized_layer_slice_[x]->Write();
            }
        }
    }
    auto fibers_Avr_LY = dirName->mkdir("fibers_Avr_LY");
    fibers_Avr_LY->cd();
    for (int j = 0; j < CHANNELS_NUMBER; ++j) {
        histograms_fiberLY_normalized_layer_[j]->Write();
    }

    auto fibers_dist_LY = dirName->mkdir("fibers_dist_LY");
    for (int ch = 0; ch < CHANNELS_NUMBER; ++ch) {
        auto x_fibers_dist_LY_N = fibers_dist_LY->mkdir(Form("x_fiber-%d_dist_LY_vs_Y", ch));
        auto y_fibers_dist_LY_N = fibers_dist_LY->mkdir(Form("y_fiber-%d_dist_LY_vs_X", ch));
        for (int y = 0; y < CUBES_SLICE_Y; ++y) {
            if (x_fiber_pixel_LY_[ch][y].size() > 10) {
                FillPixelFiberData( x_fibers_dist_LY_N, &x_fiber_pixel_LY_[ch][y], x_fiberLY_normalized_position_[ch].get(),
                x_Avr_Mean_fiberLY_normalized_position_[ch].get(), x_fiber_peak_comparison_[ch].get(), Form("fiber_x-%i_pixel_%i", ch, y), y);
            }
        }
        for (int x = 0; x < CUBES_SLICE_X; ++x) {
            if (y_fiber_pixel_LY_[ch][x].size() > 10) {
                FillPixelFiberData( y_fibers_dist_LY_N, &y_fiber_pixel_LY_[ch][x], y_fiberLY_normalized_position_[ch].get(),
                    y_Avr_Mean_fiberLY_normalized_position_[ch].get(), y_fiber_peak_comparison_[ch].get(), Form("fiber_y-%i_pixel_%i", ch, x),  x);
            }
        }
        fibers_dist_LY->cd();
        if (x_fiberLY_normalized_position_[ch]->GetEntries() > 10) {
            x_fiberLY_normalized_position_[ch]->Write();
            x_Avr_Mean_fiberLY_normalized_position_[ch]->Write();
            x_fiber_peak_comparison_[ch]->Write();
        }
        if (y_fiberLY_normalized_position_[ch]->GetEntries() > 10) {
            y_fiberLY_normalized_position_[ch]->Write();
            y_Avr_Mean_fiberLY_normalized_position_[ch]->Write();
            y_fiber_peak_comparison_[ch]->Write();
        }
    }

    dirName->cd();
    auto Gaus_layer = dirName->mkdir("Gaus_layer");
    for (int j = 0; j < 3; ++j) {
        auto gaus_layer_Summ_N = Gaus_layer->mkdir(Form("gaus_layer_Summ%d", j));
        auto gaus_layer_X_N = Gaus_layer->mkdir(Form("gaus_layer_X_%d", j));
        auto gaus_layer_Y_N = Gaus_layer->mkdir(Form("gaus_layer_Y_%d", j));
        for (int x = 0; x < CUBES_SLICE_X; ++x) {
            for (int y = 0; y < CUBES_SLICE_Y; ++y) {
                TF1* fit = new TF1("gaus","gaus", 30 , 80);
                if (!pixelsDataSummLY_[j][x][y].empty()) {
                    FillPixelData(gaus_layer_Summ_N, &pixelsDataSummLY_[j][x][y], fit, histogramsLY_gaus_layer_.at(j).get(),
                        Form("pixel_%i_%i_%i_summ", j, x, y), x, y, 200);
                }
                //----------------------------------------//----------------------------------------//
                if (!pixelsDataXLY_[j][x][y].empty()) {
                    fit = new TF1("gaus","gaus", 20 , 50);
                    FillPixelData(gaus_layer_X_N, &pixelsDataXLY_[j][x][y], fit, x_histogramsLY_gaus_layer_.at(j).get(),
                        Form("pixel_%i_%i_%i_x",j,x,y),x,y, 100);
                }
                //----------------------------------------//----------------------------------------//
                if (!pixelsDataYLY_[j][x][y].empty()) {
                    fit = new TF1("gaus","gaus", 20 , 50);
                    FillPixelData(gaus_layer_Y_N, &pixelsDataYLY_[j][x][y], fit, y_histogramsLY_gaus_layer_[j].get(),
                   Form("pixel_%i_%i_%i_y",j,x,y),x,y,100);
                }
                fit->Delete();
            }
        }
        Gaus_layer->cd();
        histogramsLY_gaus_layer_.at(j)->Write();
        x_histogramsLY_gaus_layer_.at(j)->Write();
        y_histogramsLY_gaus_layer_.at(j)->Write();
    }

    auto avr_Gaus_cube = dirName->mkdir("Avr_Gaus_cube");
    auto summ_avr_Gaus_cube = avr_Gaus_cube->mkdir("Summ_Avr_Gaus_cube");
    auto x_avr_Gaus_cube = avr_Gaus_cube->mkdir("X_Avr_Gaus_cube");
    auto y_avr_Gaus_cube = avr_Gaus_cube->mkdir("Y_Avr_Gaus_cube");
    for (int x = 0; x < CUBES_SLICE_X; ++x) {
        for (int y = 0; y < CUBES_SLICE_Y; ++y) {
            summ_avr_Gaus_cube->cd();
            TF1* fit = new TF1("gaus","gaus", 10 , 100);
            if (!pixelsAverageSummLY_[x][y].empty()) {
                FillPixelData(summ_avr_Gaus_cube, &pixelsAverageSummLY_[x][y], fit, histogramsAverageLY_gaus_layer_.get(),
                    Form("pixel_%i_%i_summ",x,y),x,y, 300);
            }
            if (!pixelsAverageXLY_[x][y].empty()) {
                fit = new TF1("gaus","gaus", 5 , 60);
                FillPixelData(x_avr_Gaus_cube, &pixelsAverageXLY_[x][y], fit, x_histogramsAverageLY_gaus_layer_.get(),
                    Form("pixel_%i_%i_x",x,y),x,y,170);
            }
            if (!pixelsAverageYLY_[x][y].empty()) {
                fit = new TF1("gaus","gaus", 5 , 60);
                FillPixelData(y_avr_Gaus_cube, &pixelsAverageYLY_[x][y], fit, y_histogramsAverageLY_gaus_layer_.get(),
                    Form("pixel_%i_%i_y",x,y),x,y,170);
            }
            fit->Delete();
        }
    }
    avr_Gaus_cube->cd();

    TCanvas tempCanvas("Gaus Mean LY X+Y", "Gaus Mean LY X+Y", 1024, 1024);
    histogramsAverageLY_gaus_layer_->Write();
    histogramsAverageLY_gaus_layer_->SetAxisRange(70, 115, "Z");
    histogramsAverageLY_gaus_layer_->SetAxisRange(-2,12,"X");
    histogramsAverageLY_gaus_layer_->SetAxisRange(-2,12,"Y");
    histogramsAverageLY_gaus_layer_->Draw("colz");
    tempCanvas.Update();
    tempCanvas.Write();
    tempCanvas.Clear();

    tempCanvas.SetName("Gaus Mean LY X");
    tempCanvas.SetTitle("Gaus Mean LY X");
    x_histogramsAverageLY_gaus_layer_->Write();
    x_histogramsAverageLY_gaus_layer_->SetAxisRange(35, 60, "Z");
    x_histogramsAverageLY_gaus_layer_->SetAxisRange(-2,12,"X");
    x_histogramsAverageLY_gaus_layer_->SetAxisRange(-2,12,"Y");
    x_histogramsAverageLY_gaus_layer_->Draw("colz");
    tempCanvas.Update();
    tempCanvas.Write();
    tempCanvas.Clear();

    tempCanvas.SetName("Gaus Mean LY Y");
    tempCanvas.SetTitle("Gaus Mean LY Y");
    y_histogramsAverageLY_gaus_layer_->Write();
    y_histogramsAverageLY_gaus_layer_->SetAxisRange(35, 60, "Z");
    y_histogramsAverageLY_gaus_layer_->SetAxisRange(-2,12,"X");
    y_histogramsAverageLY_gaus_layer_->SetAxisRange(-2,12,"Y");
    y_histogramsAverageLY_gaus_layer_->Draw("colz");
    tempCanvas.Update();
    tempCanvas.Write();
    tempCanvas.Clear();

    TSpectrum *fTSpectrum = new TSpectrum(15);
    int nfound = 0;
    TF1* fit;
    dirName->cd();
    nfound = fTSpectrum->Search(averageLY_.get(),10,"",0.5);
    double *LyPeaks = fTSpectrum->GetPositionX();
    for (uint peak = 0; peak < nfound; ++peak) {
        if (LyPeaks[peak] > 15 && LyPeaks[peak] < 75) {
            fit = new TF1("gaus","gaus", LyPeaks[peak] - 7, LyPeaks[peak] + 7);
            averageLY_->Fit(fit, "qr+");
            averageLY_->Draw();
            averageLY_->Write();
        }
    }

    nfound = fTSpectrum->Search(averageLY_X_.get(),10,"",0.5);
    LyPeaks = fTSpectrum->GetPositionX();
    for (uint peak = 0; peak < nfound; ++peak) {
        if (LyPeaks[peak] > 15 && LyPeaks[peak] < 75) {
            fit = new TF1("gaus","gaus", LyPeaks[peak] - 7, LyPeaks[peak] + 7);
            averageLY_X_->Fit(fit, "qr+");
            averageLY_X_->Draw();
            averageLY_X_->Write();
        }
    }

    nfound = fTSpectrum->Search(averageLY_Y_.get(),10,"",0.5);
    LyPeaks = fTSpectrum->GetPositionX();
    for (uint peak = 0; peak < nfound; ++peak) {
        if (LyPeaks[peak] > 15 && LyPeaks[peak] < 75) {
            fit = new TF1("gaus","gaus", LyPeaks[peak] - 7, LyPeaks[peak] + 7);
            averageLY_Y_->Fit(fit, "qr+");
            averageLY_Y_->Draw();
            averageLY_Y_->Write();
        }
    }
    averageTime_->Write();

    auto Fibers_LY_A = dirName->mkdir("Fibers_LY_A");
    Fibers_LY_A->cd();
    for (int i = 0; i < CUBES_COUNT; ++i) {
        nfound = fTSpectrum->Search(x_CubeLY_A[i].get(),10,"",0.05);
        if(nfound > 0) {
            double *xpeaks = fTSpectrum->GetPositionX();
            for (uint peak = 0; peak < nfound; ++peak) {
                if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                    fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                    x_CubeLY_A[i]->Fit(fit, "qr+");
                    break;
                }
            }
        }
        x_CubeLY_A[i]->Draw();
        x_CubeLY_A.at(i)->Write();

        nfound = fTSpectrum->Search(y_CubeLY_A[i].get(),10,"",0.05);
        if(nfound > 0) {
            double *xpeaks = fTSpectrum->GetPositionX();
            for (uint peak = 0; peak < nfound; ++peak) {
                if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                    fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                    y_CubeLY_A[i]->Fit(fit, "qr+");
                    break;
                }
            }
        }
        y_CubeLY_A[i]->Draw();
        y_CubeLY_A.at(i)->Write();
    }

    auto Fibers_LY_A_recon = dirName->mkdir("Fibers_LY_A_reconstructed");
    Fibers_LY_A_recon->cd();
    for (int i = 0; i < CUBES_COUNT; ++i) {

        nfound = fTSpectrum->Search(x_CubeLY_A_recon[i].get(),10,"",0.05);
        if(nfound > 0) {
            double *xpeaks = fTSpectrum->GetPositionX();
            for (uint peak = 0; peak < nfound; ++peak) {
                if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                    fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                    x_CubeLY_A_recon[i]->Fit(fit, "qr+");
                    break;
                }
            }
        }
        x_CubeLY_A_recon[i]->Draw();
        x_CubeLY_A_recon.at(i)->Write();

        nfound = fTSpectrum->Search(y_CubeLY_A_recon[i].get(),10,"",0.05);
        if(nfound > 0) {
            double *xpeaks = fTSpectrum->GetPositionX();
            for (uint peak = 0; peak < nfound; ++peak) {
                if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                    fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                    y_CubeLY_A_recon[i]->Fit(fit, "qr+");
                    break;
                }
            }
        }
        y_CubeLY_A_recon[i]->Draw();
        y_CubeLY_A_recon.at(i)->Write();
    }

    auto LY_1D_X_Y = dirName->mkdir(Form("LY_1D_X_Y")); // 1D
    auto LY_1D_X = LY_1D_X_Y->mkdir(Form("LY_1D_X"));
    auto LY_1D_Y = LY_1D_X_Y->mkdir(Form("LY_1D_Y"));

//    this->FillPixelAverageData_X_Y (LY_1D_X_Y, pixelsAverageSummLY_, pixelsAverageLY_X_Y_.get(),
//                                          "pixel_average_right", 200);
    this->FillPixelAverageData_X(LY_1D_X, pixelsAverageXLY_, pixelsAverageLY_X_.get(),
                                          "pixel_average_LY_X", 100);
    this->FillPixelAverageData_Y(LY_1D_Y, pixelsAverageYLY_, pixelsAverageLY_Y_.get(),
                                       "pixel_average_LY_Y", 100);


    LY_1D_X_Y->cd();
//    pixelsAverageLY_X_Y_->GetYaxis()->SetRangeUser(0, 50);
//    pixelsAverageLY_X_Y_->Write();
    pixelsAverageLY_X_->GetYaxis()->SetRangeUser(0, 50);
    pixelsAverageLY_X_->GetXaxis()->SetTitle("X or Y [mm]");
    pixelsAverageLY_X_->GetYaxis()->SetTitle("LY [p.e.]");
    pixelsAverageLY_X_->Write();

    pixelsAverageLY_Y_->GetXaxis()->SetTitle("X or Y [mm]");
    pixelsAverageLY_Y_->GetYaxis()->SetTitle("LY [p.e.]");
    pixelsAverageLY_Y_->GetYaxis()->SetRangeUser(0, 50);
    pixelsAverageLY_Y_->Write(); // 1D

    auto Fibers_LY_I = dirName->mkdir("Fibers_LY_I");
    Fibers_LY_I->cd();
    for (int i = 0; i < CUBES_COUNT; ++i) {
        nfound = fTSpectrum->Search(x_CubeLY_I[i].get(),10,"",0.05);
        if(nfound > 0) {
            double *xpeaks = fTSpectrum->GetPositionX();
            for (uint peak = 0; peak < nfound; ++peak) {
                if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                    fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                    x_CubeLY_I[i]->Fit(fit, "qr+");
                    break;
                }
            }
        }
        x_CubeLY_I[i]->Draw();
        x_CubeLY_I.at(i)->Write();

        nfound = fTSpectrum->Search(y_CubeLY_I[i].get(),10,"",0.05);
        if(nfound > 0) {
            double *xpeaks = fTSpectrum->GetPositionX();
            for (uint peak = 0; peak < nfound; ++peak) {
                if (xpeaks[peak] > 15 && xpeaks[peak] < 75) {
                    fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                    y_CubeLY_I[i]->Fit(fit, "qr+");
                    break;
                }
            }
        }
        y_CubeLY_I[i]->Draw();
        y_CubeLY_I.at(i)->Write();
    }

    auto each_Fibers_LY = dirName->mkdir("Each_Fibers_LY");
    each_Fibers_LY->cd();
    for (int ch = 0; ch < CHANNELS_NUMBER; ++ch) {
        if (fiberLY_[ch]->GetEntries() == 0 || fiberLY_[ch]->GetEntries() > 0) {
            nfound = fTSpectrum->Search(fiberLY_[ch].get(),10,"",0.05);
            if(nfound > 0) {
                double *xpeaks = fTSpectrum->GetPositionX();
                for (uint peak = 0; peak < nfound; ++peak) {
                    if (xpeaks[peak] > 15 && xpeaks[peak] < 110) {
                        if (ch <= 20){ // z волокна (с 0 по 20) идут по пучку
                            fit = new TF1("gaus","gaus", xpeaks[peak] - 15, xpeaks[peak] + 15);
                        } else {
                            fit = new TF1("gaus","gaus", xpeaks[peak] - 7, xpeaks[peak] + 7);
                        }
                        fiberLY_[ch]->Fit(fit, "qr+");
                        break;
                    }
                }
            }
            fiberLY_[ch]->Draw();
            fiberLY_[ch]->Write();
        }
    }
//     this->NormalizeCubes();
     this->NormalizeFibers();
//    averageCubeResponse_->Write();
//    averageCount_->Write();
}

void LyHistogramManager::FillPixelData(TDirectory* dir, std::vector<double>* data, TF1* fit, TH2F* hist, const char* title, const int x, const int y, const double cut) {
    dir->cd();
    const auto tempHist = new TH1F(title, title, GAUS_BINNING);
    for (const double &it: *data) {
        tempHist->Fill(it);
    }
    if (tempHist->GetEntries() > 40) {
        std::unique_ptr<TSpectrum> tSpectrum = std::make_unique<TSpectrum>(15);
//        int nFound = tSpectrum->Search(tempHist,10,"nobackground new",0.1);
//        int nFound = tSpectrum->Search(tempHist, 20, "nobackground new", 0.4);
        int nFound = tSpectrum->Search(tempHist, 2, "nobackground new", 0.05);
        if (nFound > 0) {
            double *xpeaks = tSpectrum->GetPositionX();

            int bestPeakIndex = -1;
            double maxBinContent = -1.0;

            for (uint peak = 0; peak < nFound; ++peak) {
                double x = xpeaks[peak];
                double height = tempHist->GetBinContent(tempHist->FindBin(x));

                if (height > maxBinContent) {
                    maxBinContent = height;
                    bestPeakIndex = peak;
                }
            }

            if (bestPeakIndex != -1) {
                double x_best = xpeaks[bestPeakIndex];

                if (x_best > 0 && x_best <= 3.0) { // пока спорно??
                    tempHist->Fit(fit, "qrl+", "", x_best - 0.5,x_best + 0.5); // спец. кат для очень малых световыходов
                } else if (x_best > 3.0 && x_best <= 15.0) {
                    tempHist->Fit(fit, "qrl+", "", x_best - 5.0, x_best + 5.0); // когда частица попала в волокно
                } else if (x_best > 15.0) {
                    tempHist->Fit(fit, "qrl+", "", x_best - 12.0,
                                  x_best + 12.0); // спец. кат, было 15 для большего световыхода
            }
        }
    }


        tempHist->Draw();
        tempHist->GetXaxis()->SetTitle("LY [p.e.]");
        tempHist->GetYaxis()->SetTitle("Number");
        tempHist->Write();

        double pixelLY = fit->GetParameter(1);
        double pixelSigma = fit->GetParameter(2);
        double pixelSigmaErr = fit->GetParError(2);

        double chi2 = fit->GetChisquare();
        int ndf = fit->GetNDF();

        if ((pixelLY > 0 && pixelLY < cut) && (pixelSigmaErr / pixelSigma < 0.5)){ // кат на отрезание плохого фитирования (gaus layer)
            hist->SetBinContent(x,y,pixelLY);
        }
    }
    tempHist->Delete();
}
double LyHistogramManager::FillPixelProjectionXData(TDirectory* dirName, const std::vector<double> (&pixels)[LAYERS_COUNT][CUBES_SLICE_X][CUBES_SLICE_Y], TF1* fit, TH1F* hist, const char* title, int x_, int binX_min, int binX_max, int y_, int binY_min, int binY_max, int layer, double cut) {
    double result = 0;
    const auto tempHist =  new TH1F(title, title,GAUS_BINNING);
    for (int y = binY_min; y < binY_max; ++y) {
        for (const double & it : pixels[layer][x_][y]){
            tempHist->Fill(it);
        }
    }

    if (tempHist->GetEntries() > 10) {
        std::unique_ptr<TSpectrum> tSpectrum = std::make_unique<TSpectrum>(15);
        int nFound = tSpectrum->Search(tempHist,10,"nobackground new",0.1);
        if(nFound > 0) {
            double *xpeaks = tSpectrum->GetPositionX();
            for (uint peak = 0; peak < nFound; ++peak) {
                if (xpeaks[peak] > 15) {
                    tempHist->Fit(fit,"qr+", "",xpeaks[peak] - 15, xpeaks[peak] + 15);
                    break;
                }
            }
        }

        double pixelLY = fit->GetParameter(1);
        if (pixelLY > 0) {
            hist->SetBinContent(x_ - binX_min, pixelLY);
            result = pixelLY;
        }
    }
    dirName->cd();
    // tempHist->Write();

    delete tempHist;
    return  result;
}

double LyHistogramManager::FillPixelProjectionYData(TDirectory* dirName, const std::vector<double> (&pixels)[LAYERS_COUNT][CUBES_SLICE_X][CUBES_SLICE_Y], TF1* fit, TH1F* hist, const char* title, int x_, int binX_min, int binX_max, int y_, int binY_min, int binY_max, int layer, double cut) {
    double result = 0;
    const auto tempHist =  new TH1F(title, title,GAUS_BINNING);
    for (int x = binX_min; x < binX_max; ++x) {
        for (const double & it : pixels[layer][x][y_]){
            tempHist->Fill(it);
        }
    }

    if (tempHist->GetEntries() > 10) {
        std::unique_ptr<TSpectrum> tSpectrum = std::make_unique<TSpectrum>(15);
        int nFound = tSpectrum->Search(tempHist,10,"nobackground new",0.1);
        if(nFound > 0) {
            double *xpeaks = tSpectrum->GetPositionX();
            for (uint peak = 0; peak < nFound; ++peak) {
                if (xpeaks[peak] > 15) {
                    tempHist->Fit(fit,"qr+", "",xpeaks[peak] - 15, xpeaks[peak] + 15);
                    break;
                }
            }
        }

        double pixelLY = fit->GetParameter(1);
        if (pixelLY > 0 ) {
            hist->SetBinContent(y_ - binY_min, pixelLY);
            result = pixelLY;
        }
    }
    dirName->cd();
    // tempHist->Write();

    delete tempHist;
    return  result;
}

double LyHistogramManager::FillPixelData(std::vector<double>* data, TF1* fit, TH2F* hist, TH2F* histError, TH2F* histSigma, const char* title, int x, int y, double cut) {
    double result = 0;
    const auto tempHist =  new TH1F(title, title, GAUS_BINNING);
    for (const double & it : *data) {
        tempHist->Fill(it);
    }
    if (tempHist->GetEntries() > 10) {
        std::unique_ptr<TSpectrum> tSpectrum = std::make_unique<TSpectrum>(15);
//        int nFound = tSpectrum->Search(tempHist,10,"nobackground new",0.1);
//        int nFound = tSpectrum->Search(tempHist, 20, "nobackground new", 0.4);
        int nFound = tSpectrum->Search(tempHist, 2, "nobackground new", 0.05);
        if(nFound > 0) {
            double *xpeaks = tSpectrum->GetPositionX();
            for (uint peak = 0; peak < nFound; ++peak) {
                if (xpeaks[peak] > 15) {
                    tempHist->Fit(fit,"qr+", "",xpeaks[peak] - 15, xpeaks[peak] + 15);
                    break;
                }
            }
        }

        double pixelLY = fit->GetParameter(1);
        double pixelLYError = fit->GetParError(1);
        double pixelSigma = fit->GetParameter(2);
        double pixelSigmaErr = fit->GetParError(2);

        double chi2 = fit->GetChisquare();
        int ndf = fit->GetNDF();

        if ((pixelLY > 0 && pixelLY < cut && pixelLYError < 20) && ndf != 0 && (pixelSigmaErr / pixelSigma < 0.5)) { // кат на отрезание плохого фитирования
            hist->SetBinContent(x,y,pixelLY);
            histError->SetBinContent(x,y,pixelLYError);
            histSigma->SetBinContent(x,y,pixelSigma);
            result = pixelLY;
        }
    }
    tempHist->Delete();
    return result;
}

void LyHistogramManager::FillPixelFiberData(TDirectory* dir, std::vector<double>* data, TH1F* hist, TH1F* histAvr, TH1F* histCompr, const char* title, int position) {
    dir->cd();
    TH1F *tempHist = new TH1F(title, title, GAUS_FIBER_BINNING);
    unsigned int peak1_count = 0;
    unsigned int peak2_count = 0;

    for (const double &it: *data) {
        tempHist->Fill(it);
        it < 18 ? ++peak1_count : ++peak2_count;
    }

    if (tempHist->GetEntries() > 500) { // спец. кат, раньше был 15000, далее был 5000

        std::unique_ptr<TSpectrum> tSpectrum = std::make_unique<TSpectrum>(15);
        int nFound = tSpectrum->Search(tempHist, 10, "nobackground new", 0.1);
        TF1 *fit = new TF1("gaus", "gaus");
        if (nFound > 0) {
            double *xpeaks = tSpectrum->GetPositionX();
            for (uint peak = 0; peak < nFound; ++peak) {
                if (xpeaks[peak] >= 18 && xpeaks[peak] < 70) {
                    fit = new TF1("gaus", "gaus");
                    tempHist->Fit(fit, "qr+", "", xpeaks[peak] - 15, xpeaks[peak] + 15);
                    break;
                }
                if (xpeaks[peak] > 0 && xpeaks[peak] < 18) {
                    fit = new TF1("landau", "landau");
                    tempHist->Fit(fit, "qr+", "", xpeaks[peak] - 5, xpeaks[peak] + 5);
                    break;
                }
            }
        }
        tempHist->Draw();
        tempHist->GetXaxis()->SetTitle("LY [p.e.]");
        tempHist->GetYaxis()->SetTitle("Number");
        tempHist->Write();
        // std::cout << "total: " << data->size() << " peak_1: " << peak1_count << " peak_2: " << peak2_count << std::endl;
        if (tempHist->GetEntries() > 400)
            histCompr->SetBinContent(position, static_cast<double>(peak1_count) / static_cast<double>(peak2_count));

        double pixelLY = fit->GetParameter(1);
        if (pixelLY > 0 && pixelLY < 80) {
            hist->SetBinContent(position, pixelLY);
        }

        histAvr->SetBinContent(position, tempHist->GetMean());
        tempHist->Delete();
    }
}


void LyHistogramManager::FillPixelAverageData_X(TDirectory* dir, const std::vector<double>(&data)[CUBES_SLICE_X][CUBES_SLICE_Y], TH1F* hist, const std::string& title, double cut) { // 1D
    dir->cd();
    std::unique_ptr<TH1F> tempHist[CUBES_SLICE_X];
    for (int x = 0; x < CUBES_SLICE_X; ++x) {
        std::string name = title + '_' + std::to_string(x);
        tempHist[x] =  std::make_unique<TH1F>(name.c_str(), name.c_str(), LY_BINNING_2);
    }
    for (int x = 0; x < CUBES_SLICE_X; ++x) {
        for (int y = 0; y < CUBES_SLICE_Y; ++y) {
            for (const auto it : data[x][y]) {
                tempHist[x]->Fill(it);
            }
        }

        if (tempHist[x]->GetEntries() > 500) { // спец. кат
            tempHist[x]->GetXaxis()->SetTitle("LY [p.e.]");
            tempHist[x]->GetYaxis()->SetTitle("N");
            TF1* gausFit = new TF1("gaus", "gaus", 15, 35);
            tempHist[x]->Fit(gausFit, "QR");
            tempHist[x]->Draw();
            tempHist[x]->Write();
//            double pixel_mean = tempHist[x]->GetMean();
            double pixel_mean = gausFit->GetParameter(1);
            if ((pixel_mean == 0 ||  pixel_mean > 0) && pixel_mean < cut) {
                hist->SetBinContent(x, pixel_mean);
            }
        }
    }
}

void LyHistogramManager::FillPixelAverageData_Y(TDirectory* dir, const std::vector<double>(&data)[CUBES_SLICE_X][CUBES_SLICE_Y], TH1F* hist, const std::string& title, double cut) {
    dir->cd();
    std::unique_ptr<TH1F> tempHist[CUBES_SLICE_Y];
    for (int y = 0; y < CUBES_SLICE_Y; ++y) {
        std::string name = title + '_' + std::to_string(y);
        tempHist[y] =  std::make_unique<TH1F>(name.c_str(), name.c_str(), LY_BINNING_2);
    }
    for (int y = 0; y < CUBES_SLICE_Y; ++y) {
        for (int x = 0; x < CUBES_SLICE_X; ++x) {
            for (const auto it : data[x][y]) {
                tempHist[y]->Fill(it);
            }
        }
        if ( tempHist[y]->GetEntries() > 500) { // спец. кат
            tempHist[y]->GetXaxis()->SetTitle("LY [p.e.]");
            tempHist[y]->GetYaxis()->SetTitle("N");
            TF1* gausFit = new TF1("gaus", "gaus", 15, 35);
            tempHist[y]->Fit(gausFit, "QR+");
            tempHist[y]->Draw();
            tempHist[y]->Write();
//            double pixel_mean = tempHist[y]->GetMean();
//            double pixel_max = tempHist[y]->GetMaximum();
            double pixel_mean = gausFit->GetParameter(1);
            if (pixel_mean > 0 && pixel_mean < cut) {
                hist->SetBinContent(y, pixel_mean);
            }
        }
    }
} // 1D

XTalkHistogramManager::XTalkHistogramManager(unsigned int count) : HistogramManager(count){
    x_histogramsXTalkN_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_XTalk_N_" + std::to_string(i);
        std::string title = "x_XTalk_N_" + std::to_string(i);
        x_histogramsXTalkN_.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), CUBES_BINNING_X));
    }

    y_histogramsXTalkN_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_XTalk_N_" + std::to_string(i);
        std::string title = "y_XTalk_N_" + std::to_string(i);
        y_histogramsXTalkN_.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), CUBES_BINNING_Y));
    }

    x_histogramsXTalk_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_XTalk_" + std::to_string(i);
        std::string title = "x_XTalk_" + std::to_string(i);
        x_histogramsXTalk_.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), CUBES_BINNING_X));
    }

    y_histogramsXTalk_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_XTalk_" + std::to_string(i);
        std::string title = "y_XTalk_" + std::to_string(i);
        y_histogramsXTalk_.emplace_back(std::make_unique<TH1F>(name.c_str(), title.c_str(), CUBES_BINNING_Y));
    }

    histogramsXTalk_2D_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "XTalk_2D_" + std::to_string(i);
        std::string title = "XTalk_2D_" + std::to_string(i);
        histogramsXTalk_2D_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    x_histogramsXTalk_2D_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_XTalk_2D_" + std::to_string(i);
        std::string title = "x_XTalk_2D_" + std::to_string(i);
        x_histogramsXTalk_2D_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    y_histogramsXTalk_2D_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_XTalk_2D_" + std::to_string(i);
        std::string title = "y_XTalk_2D_" + std::to_string(i);
        y_histogramsXTalk_2D_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    x_histogramsLY_gaus_cube_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "x_XTalk_2D_Gaus" + std::to_string(i);
        std::string title = "x_XTalk_2D_Gaus" + std::to_string(i);
        x_histogramsLY_gaus_cube_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    y_histogramsLY_gaus_cube_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        std::string name = "y_XTalk_2D_Gaus" + std::to_string(i);
        std::string title = "y_XTalk_2D_Gaus" + std::to_string(i);
        y_histogramsLY_gaus_cube_.emplace_back(std::make_unique<TH2F>(name.c_str(), title.c_str(), CUBES_BINNING));
    }

    averageXTalk_ = std::make_unique<TH1F>("Average_XTalk", "Average_XTalk", 200, 0 , 1);
    averageXTalk_x_ = std::make_unique<TH1F>("Average_XTalk_X", "Average_XTalk_X", 200, 0 , 1);
    averageXTalk_y_ = std::make_unique<TH1F>("Average_XTalk_Y", "Average_XTalk_Y", 200, 0 , 1);

    pixelsDataXLY_cube_.resize(CUBES_COUNT);
    pixelsDataYLY_cube_.resize(CUBES_COUNT);
    for (int cube = 0; cube < CUBES_COUNT; ++cube) {
        pixelsDataXLY_cube_[cube].resize(CUBES_SLICE_X);
        pixelsDataYLY_cube_[cube].resize(CUBES_SLICE_X);
        for (int x = 0; x < CUBES_SLICE_X; ++x) {
            pixelsDataXLY_cube_[cube][x].resize(CUBES_SLICE_Y);
            pixelsDataYLY_cube_[cube][x].resize(CUBES_SLICE_Y);
        }
    }
}

void XTalkHistogramManager::SetTitlesHistXTalk_N(const std::string& titleX, const std::string& titleY) {
    for (auto &a: x_histogramsXTalkN_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
    }

    for (auto &a: y_histogramsXTalkN_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
    }
}

void XTalkHistogramManager::SetTitlesHistXTalk_1D(const std::string& titleX, const std::string& titleY){
    for (auto& a : x_histogramsXTalk_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
    }

    for (auto& a : y_histogramsXTalk_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
    }
}

void XTalkHistogramManager::SetTitlesHistXTalk_1D_all(const std::string& titleX, const std::string& titleY){
    if (averageXTalk_) {
        averageXTalk_->GetXaxis()->SetTitle(titleX.c_str());
        averageXTalk_->GetYaxis()->SetTitle(titleY.c_str());
    }

    if (averageXTalk_x_) {
        averageXTalk_x_->GetXaxis()->SetTitle(titleX.c_str());
        averageXTalk_x_->GetYaxis()->SetTitle(titleY.c_str());
    }

    if (averageXTalk_y_) {
        averageXTalk_y_->GetXaxis()->SetTitle(titleX.c_str());
        averageXTalk_y_->GetYaxis()->SetTitle(titleY.c_str());
    }

}

void XTalkHistogramManager::SetTitlesHistXTalk_2D(const std::string& titleX, const std::string& titleY, const std::string& titleZ){
    for (auto& a : histogramsXTalk_2D_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle((titleZ).c_str());
    }

    for (auto& a : x_histogramsXTalk_2D_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle((titleZ).c_str());
    }

    for (auto& a : y_histogramsXTalk_2D_) {
        a->GetXaxis()->SetTitle((titleX).c_str());
        a->GetYaxis()->SetTitle((titleY).c_str());
        a->GetZaxis()->SetTitle((titleZ).c_str());
    }
}

void XTalkHistogramManager::SetBinContentXTalk_x(unsigned int index, double x, double y, double x_value, double y_value){
    //---------1D----------//
    int x_binXly = x_histogramsXTalk_[index]->GetXaxis()->FindBin(x);
    auto x_currentValue_norm = x_histogramsXTalk_[index]->GetBinContent(x_binXly);

    int binX = x_histogramsXTalkN_[index]->GetXaxis()->FindBin(x);
    auto currentCountX = x_histogramsXTalkN_[index]->GetBinContent(binX);

    x_histogramsXTalkN_[index]->Fill(x);

    double x_newAvr_norm = (x_currentValue_norm * currentCountX + x_value)/(currentCountX + 1);
    x_histogramsXTalk_[index]->SetBinContent(x_binXly, x_newAvr_norm);

    //---------2D----------//
    int x_binXly_2D = x_histogramsXTalk_2D_[index]->GetXaxis()->FindBin(x);
    int y_binXly_2D = x_histogramsXTalk_2D_[index]->GetYaxis()->FindBin(y);
    auto x_currentValue_norm_2D = x_histogramsXTalk_2D_[index]->GetBinContent(x_binXly_2D, y_binXly_2D);

    int binX_2D = histogramsXTalk_2D_[index]->GetXaxis()->FindBin(x);
    int binY_2D = histogramsXTalk_2D_[index]->GetYaxis()->FindBin(y);
    auto currentCountX_2D = histogramsXTalk_2D_[index]->GetBinContent(binX_2D, binY_2D);

    histogramsXTalk_2D_[index]->Fill(x, y);

    double x_newAvr_norm_2D = (x_currentValue_norm_2D * currentCountX_2D + x_value)/(currentCountX_2D + 1);
    x_histogramsXTalk_2D_[index]->SetBinContent(x_binXly_2D, y_binXly_2D, x_newAvr_norm_2D);

    averageXTalk_->Fill(x_value);
    averageXTalk_x_->Fill(x_value);

    pixelsDataXLY_cube_[index][x_binXly_2D][y_binXly_2D].push_back(x_value);
}

void XTalkHistogramManager::SetBinContentXTalk_y(unsigned int index, double x, double y, double x_value, double y_value){
    //---------1D----------//
    int y_binXly = y_histogramsXTalk_[index]->GetXaxis()->FindBin(y);
    auto y_currentValue_norm = y_histogramsXTalk_[index]->GetBinContent(y_binXly);


    int binY = y_histogramsXTalkN_[index]->GetXaxis()->FindBin(y);
    auto currentCountY = y_histogramsXTalkN_[index]->GetBinContent(binY);

    y_histogramsXTalkN_[index]->Fill(y);

    double y_newAvr_norm = (y_currentValue_norm * currentCountY + y_value)/(currentCountY + 1);
    y_histogramsXTalk_[index]->SetBinContent(y_binXly, y_newAvr_norm);

    //---------2D----------//
    int x_binXly_2D = y_histogramsXTalk_2D_[index]->GetXaxis()->FindBin(x);
    int y_binXly_2D = y_histogramsXTalk_2D_[index]->GetYaxis()->FindBin(y);
    auto y_currentValue_norm_2D = y_histogramsXTalk_2D_[index]->GetBinContent(x_binXly_2D, y_binXly_2D);

    int binX_2D = histogramsXTalk_2D_[index]->GetXaxis()->FindBin(x);
    int binY_2D = histogramsXTalk_2D_[index]->GetYaxis()->FindBin(y);
    auto currentCountY_2D = histogramsXTalk_2D_[index]->GetBinContent(binX_2D, binY_2D);

    histogramsXTalk_2D_[index]->Fill(x, y);

    double y_newAvr_norm_2D = (y_currentValue_norm_2D * currentCountY_2D + y_value)/(currentCountY_2D + 1);
    y_histogramsXTalk_2D_[index]->SetBinContent(x_binXly_2D, y_binXly_2D, y_newAvr_norm_2D);

    averageXTalk_->Fill(y_value);
    averageXTalk_y_->Fill(y_value);

    pixelsDataYLY_cube_[index][x_binXly_2D][y_binXly_2D].push_back(y_value);
}

void XTalkHistogramManager::WriteHistograms(TDirectory* dirName){
    auto avrXCount_x = dirName->mkdir("avrXCount_x");
    avrXCount_x->cd();
    for (const auto& a : x_histogramsXTalkN_) {
        a->Write();
    }

    auto avrXCount_y = dirName->mkdir("avrXCount_y");
    avrXCount_y->cd();
    for (const auto& a : y_histogramsXTalkN_) {
        a->Write();
    }

    auto avrXTalk_x = dirName->mkdir("avrXTalk_1D_x");
    avrXTalk_x->cd();
    for (const auto& a : x_histogramsXTalk_) {
        a->Write();
    }

    auto avrXTalk_y = dirName->mkdir("avrXTalk_1D_y");
    avrXTalk_y->cd();
    for (const auto& a : y_histogramsXTalk_) {
        a->Write();
    }

    auto avrXCount_2D = dirName->mkdir("avrXCount_2D");
    avrXCount_2D->cd();
    for (const auto& a : histogramsXTalk_2D_) {
        a->Write();
    }

    auto avrXTalk_2D_x = dirName->mkdir("avrXTalk_2D_x");
    avrXTalk_2D_x->cd();
    for (const auto& a : x_histogramsXTalk_2D_) {
        a->Write();
    }

    auto avrXTalk_2D_y = dirName->mkdir("avrXTalk_2D_y");
    avrXTalk_2D_y->cd();
    for (const auto& a : y_histogramsXTalk_2D_) {
        a->Write();
    }

//    TF1 *landauFit = new TF1("landauFit", "landau", 0, 1);

    auto average_XTalk_ = dirName->mkdir("average_XTalk_");
    average_XTalk_->cd();

//    averageXTalk_->Fit("landau", "r", " ", 0.02, 0.5);
    averageXTalk_->Draw();
    averageXTalk_->Write();

//    averageXTalk_x_->Fit("landau", "r", " ", 0.03, 0.5);
    averageXTalk_x_->Draw();
    averageXTalk_x_->Write();

//    averageXTalk_y_->Fit("landau", "r", " ", 0.02, 0.5);
    averageXTalk_y_->Draw();
    averageXTalk_y_->Write();

    auto Gaus_layer = dirName->mkdir("Gaus_cubes");
    for (int cube = 0; cube < CUBES_COUNT; ++cube) {
        auto gaus_layer_X_N = Gaus_layer->mkdir(Form("gaus_cube_%d_X", cube));
        auto gaus_layer_Y_N = Gaus_layer->mkdir(Form("gaus_cybe_%d_Y", cube));
        for (int x = 0; x < CUBES_SLICE_X; ++x) {
            for (int y = 0; y < CUBES_SLICE_Y; ++y) {
//                TF1* fit = new TF1("landau","landau",0, 0.3); // переход от Ландау к Mean
                double pixel_mean = 0.0;
                if (pixelsDataXLY_cube_[cube][x][y].size() > 0) {
                    gaus_layer_X_N->cd();
                    auto tempHist_X =  new TH1F(Form("pixel_cube_%i_x-%i_y-%i_X",cube,x,y),Form("pixel_cube_%i_x-%i_y-%i_X", cube, x, y), XTALK_BINNING);
                    for (auto it = pixelsDataXLY_cube_[cube][x][y].begin(); it != pixelsDataXLY_cube_[cube][x][y].end(); ++it) {
                        tempHist_X->Fill(*it);
                    }
                    if (tempHist_X->GetEntries()>0) {
//                        tempHist_X->Fit(fit, "qr+");
                        pixel_mean = tempHist_X->GetMean();
                        // double pixel_mean_err = tempHist_X->GetMeanError(); // ошибка среднего
                        tempHist_X->Draw();
                        tempHist_X->Write();
                    }
//                    pixelLY = fit->GetParameter(1);
                    if (pixel_mean >= 0 && pixel_mean < 0.5) {
                        x_histogramsLY_gaus_cube_.at(cube)->SetBinContent(x,y,pixel_mean);
                    }
                    tempHist_X->Delete();
                }
                //----------------------------------------//----------------------------------------//
                if (pixelsDataYLY_cube_[cube][x][y].size() > 0) {
                    gaus_layer_Y_N->cd();
                    auto tempHist_Y =  new TH1F(Form("pixel_cube_%i_x-%i_y-%i_Y",cube, x, y),Form("pixel_cube_%i_x-%i_y-%i_Y", cube, x, y), XTALK_BINNING);
                    for (auto it = pixelsDataYLY_cube_[cube][x][y].begin(); it != pixelsDataYLY_cube_[cube][x][y].end(); ++it) {
                        tempHist_Y->Fill(*it);
                    }
                    if (tempHist_Y->GetEntries()>0) {
//                        tempHist_Y->Fit(fit, "qr+");
                        tempHist_Y->Draw();
                        tempHist_Y->Write();
                    }
                    pixel_mean = tempHist_Y->GetMean();
                    // double pixel_mean_err = tempHist_X->GetMeanError(); // ошибка среднего
                    if (pixel_mean >= 0 && pixel_mean < 0.5) {
                        y_histogramsLY_gaus_cube_.at(cube)->SetBinContent(x, y, pixel_mean);
                    }
                    tempHist_Y->Delete();
                }
//                fit->Delete();
            }
        }
        Gaus_layer->cd();
        x_histogramsLY_gaus_cube_.at(cube)->Write();
        y_histogramsLY_gaus_cube_.at(cube)->Write();
    }
}

// новое, касательно 4 сторон
XTalk4sideHistogramManager::XTalk4sideHistogramManager(unsigned int count) : HistogramManager(count){
    histogramsXTalk_2D_up_ = std::make_unique<TH2F>("Average_XTalk_2D_Up", "Average_XTalk_2D_Up", CUBES_BINNING);
    histogramsXTalk_2D_down_ = std::make_unique<TH2F>("Average_XTalk_2D_Down", "Average_XTalk_2D_Down", CUBES_BINNING);
    histogramsXTalk_2D_right_ = std::make_unique<TH2F>("Average_XTalk_2D_Right", "Average_XTalk_2D_Right", CUBES_BINNING);
    histogramsXTalk_2D_left_ = std::make_unique<TH2F>("Average_XTalk_2D_Left", "Average_XTalk_2D_Left", CUBES_BINNING);
    histogramsXTalk_2D_horiz_ = std::make_unique<TH2F>("Average_XTalk_2D_Horizontal", "Average_XTalk_2D_Horizontal", CUBES_BINNING);
    histogramsXTalk_2D_vert_ = std::make_unique<TH2F>("Average_XTalk_2D_Vertical", "Average_XTalk_2D_Vertical", CUBES_BINNING);

    averageXTalk_up_ = std::make_unique<TH1F>("Average_XTalk_Up", "Average_XTalk_Up", 200, 0, 1);
    averageXTalk_down_ = std::make_unique<TH1F>("Average_XTalk_Down", "Average_XTalk_Down", 200, 0, 1);
    averageXTalk_right_ = std::make_unique<TH1F>("Average_XTalk_Right", "Average_XTalk_Right", 200, 0, 1);
    averageXTalk_left_ = std::make_unique<TH1F>("Average_XTalk_Left", "Average_XTalk_Left", 200, 0, 1);
    averageXTalk_horiz_ = std::make_unique<TH1F>("Average_XTalk_Horizontal", "Average_XTalk_Horizontal", 200, 0, 1);
    averageXTalk_vert_ = std::make_unique<TH1F>("Average_XTalk_Vertical", "Average_XTalk_Vertical", 200, 0, 1);

    dist_LY_XTalk_up_ = std::make_unique<TH1F>("LY_XTalk_Up", "LY_XTalk_Up", 2 * 150, 0, 150);
    dist_LY_XTalk_down_ = std::make_unique<TH1F>("LY_XTalk_Down", "LY_XTalk_Down", 2 * 150, 0, 150);
    dist_LY_XTalk_right_ = std::make_unique<TH1F>("LY_XTalk_Right", "LY_XTalk_Right", 2 * 150, 0, 150);
    dist_LY_XTalk_left_ = std::make_unique<TH1F>("LY_XTalk_Left", "LY_XTalk_Left", 2 * 150, 0, 150);
    dist_LY_XTalk_horiz_ = std::make_unique<TH1F>("LY_XTalk_Horizontal", "LY_XTalk_Horizontal", 2 * 150, 0, 150);
    dist_LY_XTalk_vert_ = std::make_unique<TH1F>("LY_XTalk_Vertical", "LY_XTalk_Vertical", 2 * 150, 0, 150);

    dist_LY_up_ = std::make_unique<TH1F>("LY_Up", "LY_Up", 5000, 0, 200);
    dist_LY_down_ = std::make_unique<TH1F>("LY_Down", "LY_Down", 5000, 0, 200);
    dist_LY_right_ = std::make_unique<TH1F>("LY_Right", "LY_Right", 5000, 0, 200);
    dist_LY_left_ = std::make_unique<TH1F>("LY_Left", "LY_Left", 5000, 0, 200);
    dist_LY_horiz_ = std::make_unique<TH1F>("LY_Horizontal", "LY_Horizontal", 5000, 0, 200);
    dist_LY_vert_ = std::make_unique<TH1F>("LY_Vertical", "LY_Vertical", 5000, 0, 200);

    pixelsAverageXTalk_up_ = std::make_unique<TH1F>("Pixel_Average_XTalk_Up", "Pixel_Average_XTalk_Up", CUBES_BINNING_Y);
    pixelsAverageXTalk_down_ = std::make_unique<TH1F>("Pixel_Average_XTalk_Down", "Pixel_Average_XTalk_Down", CUBES_BINNING_Y);
    pixelsAverageXTalk_right_ = std::make_unique<TH1F>("Pixel_Average_XTalk_Right", "Pixel_Average_XTalk_Right", CUBES_BINNING_X);
    pixelsAverageXTalk_left_ = std::make_unique<TH1F>("Pixel_Average_XTalk_Left", "Pixel_Average_XTalk_Left", CUBES_BINNING_X);
    pixelsAverageXTalk_horiz_ = std::make_unique<TH1F>("Pixel_Average_XTalk_Horizontal", "Pixel_Average_XTalk_Horizontal", CUBES_BINNING_X);
    pixelsAverageXTalk_vert_ = std::make_unique<TH1F>("Pixel_Average_XTalk_Vertical", "Pixel_Average_XTalk_Vertical", CUBES_BINNING_Y);
}

void XTalk4sideHistogramManager::SetTitlesHistXTalk4side_1D(const std::string& titleX, const std::string& titleY){
    std::array<std::unique_ptr<TH1F>*, 6> histograms = {&averageXTalk_up_, &averageXTalk_down_, &averageXTalk_right_, &averageXTalk_left_, &averageXTalk_horiz_, &averageXTalk_vert_};

    for (auto* histPtr : histograms) {
        if (histPtr && *histPtr) {
            (*histPtr)->GetXaxis()->SetTitle(titleX.c_str());
            (*histPtr)->GetYaxis()->SetTitle(titleY.c_str());
        }
    }
}

void XTalk4sideHistogramManager::SetTitlesHistXTalk4side_1D_2(const std::string& titleX, const std::string& titleY) {
    std::array<std::unique_ptr<TH1F>*, 4> histograms = {&pixelsAverageXTalk_right_, &pixelsAverageXTalk_left_, &pixelsAverageXTalk_up_, &pixelsAverageXTalk_down_};

    for (auto* histPtr : histograms) {
        if (histPtr && *histPtr) {
            (*histPtr)->GetXaxis()->SetTitle(titleX.c_str());
            (*histPtr)->GetYaxis()->SetTitle(titleY.c_str());
        }
    }
}

void XTalk4sideHistogramManager::SetTitlesHistXTalk4side_1D_3(const std::string& titleX, const std::string& titleY) {
    std::array<std::unique_ptr<TH1F>*, 11> histograms = {&dist_LY_XTalk_up_, &dist_LY_XTalk_down_, &dist_LY_XTalk_right_, &dist_LY_XTalk_left_, &dist_LY_XTalk_horiz_, &dist_LY_XTalk_vert_, &dist_LY_up_, &dist_LY_down_, &dist_LY_right_, &dist_LY_left_, &dist_LY_horiz_};

    for (auto* histPtr : histograms) {
        if (histPtr && *histPtr) {
            (*histPtr)->GetXaxis()->SetTitle(titleX.c_str());
            (*histPtr)->GetYaxis()->SetTitle(titleY.c_str());
        }
    }
}

void XTalk4sideHistogramManager::SetTitlesHistXTalk4side_2D(const std::string &titleX, const std::string &titleY, const std::string &titleZ) {
    std::array<std::unique_ptr<TH2F>*, 6> histograms = {&histogramsXTalk_2D_up_, &histogramsXTalk_2D_down_, &histogramsXTalk_2D_right_, &histogramsXTalk_2D_left_, &histogramsXTalk_2D_horiz_, &histogramsXTalk_2D_vert_};

    for (auto* histPtr : histograms) {
        if (histPtr && *histPtr) {
            (*histPtr)->GetXaxis()->SetTitle(titleX.c_str());
            (*histPtr)->GetYaxis()->SetTitle(titleY.c_str());
            (*histPtr)->GetZaxis()->SetTitle(titleZ.c_str());
        }
    }
}

void XTalk4sideHistogramManager::SetBinContentXTalk4side_up(double x, double y, double value_XTalk, double value_LY, double value_LY_XTalk){
    //---------2D----------//
    int x_binXly_2D = histogramsXTalk_2D_up_->GetXaxis()->FindBin(x);
    int y_binXly_2D = histogramsXTalk_2D_up_->GetYaxis()->FindBin(y);
    auto x_currentValue_norm_2D = histogramsXTalk_2D_up_->GetBinContent(x_binXly_2D, y_binXly_2D);

    int binX_2D = histogramsXTalk_2D_up_->GetXaxis()->FindBin(x);
    int binY_2D = histogramsXTalk_2D_up_->GetYaxis()->FindBin(y);
    auto currentCountX_2D = histogramsXTalk_2D_up_->GetBinContent(binX_2D, binY_2D);

//    histogramsXTalk_2D_up_->Fill(x, y);

    double x_newAvr_norm_2D = (x_currentValue_norm_2D * currentCountX_2D + value_XTalk)/(currentCountX_2D + 1);
//    histogramsXTalk_2D_up_->SetBinContent(x_binXly_2D, y_binXly_2D, x_newAvr_norm_2D);

    averageXTalk_up_->Fill(value_XTalk);
    dist_LY_XTalk_up_->Fill(value_LY_XTalk);
    dist_LY_up_->Fill(value_LY);
    pixelsData_up_LY_cube_[x_binXly_2D][y_binXly_2D].push_back(value_XTalk);
}

void XTalk4sideHistogramManager::SetBinContentXTalk4side_down(double x, double y,double value_XTalk, double value_LY, double value_LY_XTalk){
    //---------2D----------//
    int x_binXly_2D = histogramsXTalk_2D_down_->GetXaxis()->FindBin(x);
    int y_binXly_2D = histogramsXTalk_2D_down_->GetYaxis()->FindBin(y);
    auto x_currentValue_norm_2D = histogramsXTalk_2D_down_->GetBinContent(x_binXly_2D, y_binXly_2D);

    int binX_2D = histogramsXTalk_2D_down_->GetXaxis()->FindBin(x);
    int binY_2D = histogramsXTalk_2D_down_->GetYaxis()->FindBin(y);
    auto currentCountX_2D = histogramsXTalk_2D_down_->GetBinContent(binX_2D, binY_2D);

//    histogramsXTalk_2D_down_->Fill(x, y);

    double x_newAvr_norm_2D = (x_currentValue_norm_2D * currentCountX_2D + value_XTalk)/(currentCountX_2D + 1);
//    histogramsXTalk_2D_down_->SetBinContent(x_binXly_2D, y_binXly_2D, x_newAvr_norm_2D);

    averageXTalk_down_->Fill(value_XTalk);
    dist_LY_XTalk_down_->Fill(value_LY_XTalk);
    dist_LY_down_->Fill(value_LY);
    pixelsData_down_LY_cube_[x_binXly_2D][y_binXly_2D].push_back(value_XTalk);
}

void XTalk4sideHistogramManager::SetBinContentXTalk4side_right(double x, double y, double value_XTalk, double value_LY, double value_LY_XTalk){
    //---------2D----------//
    int x_binXly_2D = histogramsXTalk_2D_right_->GetXaxis()->FindBin(x);
    int y_binXly_2D = histogramsXTalk_2D_right_->GetYaxis()->FindBin(y);
    auto x_currentValue_norm_2D = histogramsXTalk_2D_right_->GetBinContent(x_binXly_2D, y_binXly_2D);

    int binX_2D = histogramsXTalk_2D_right_->GetXaxis()->FindBin(x);
    int binY_2D = histogramsXTalk_2D_right_->GetYaxis()->FindBin(y);
    auto currentCountX_2D = histogramsXTalk_2D_right_->GetBinContent(binX_2D, binY_2D);

//    histogramsXTalk_2D_right_->Fill(x, y);

    double x_newAvr_norm_2D = (x_currentValue_norm_2D * currentCountX_2D + value_XTalk)/(currentCountX_2D + 1);
//    histogramsXTalk_2D_right_->SetBinContent(x_binXly_2D, y_binXly_2D, x_newAvr_norm_2D);

    averageXTalk_right_->Fill(value_XTalk);
    dist_LY_XTalk_right_->Fill(value_LY_XTalk);
    dist_LY_right_->Fill(value_LY);
    pixelsData_right_LY_cube_[x_binXly_2D][y_binXly_2D].push_back(value_XTalk);
}

void XTalk4sideHistogramManager::SetBinContentXTalk4side_left(double x, double y, double value_XTalk, double value_LY, double value_LY_XTalk){
    //---------2D----------//
    int x_binXly_2D = histogramsXTalk_2D_left_->GetXaxis()->FindBin(x);
    int y_binXly_2D = histogramsXTalk_2D_left_->GetYaxis()->FindBin(y);
    auto x_currentValue_norm_2D = histogramsXTalk_2D_left_->GetBinContent(x_binXly_2D, y_binXly_2D);

    int binX_2D = histogramsXTalk_2D_left_->GetXaxis()->FindBin(x);
    int binY_2D = histogramsXTalk_2D_left_->GetYaxis()->FindBin(y);
    auto currentCountX_2D = histogramsXTalk_2D_left_->GetBinContent(binX_2D, binY_2D);

//    histogramsXTalk_2D_left_->Fill(x, y);

    double x_newAvr_norm_2D = (x_currentValue_norm_2D * currentCountX_2D + value_XTalk)/(currentCountX_2D + 1);
//    histogramsXTalk_2D_left_->SetBinContent(x_binXly_2D, y_binXly_2D, x_newAvr_norm_2D);

    averageXTalk_left_->Fill(value_XTalk);
    dist_LY_XTalk_left_->Fill(value_LY_XTalk);
    dist_LY_left_->Fill(value_LY);
    pixelsData_left_LY_cube_[x_binXly_2D][y_binXly_2D].push_back(value_XTalk);
}

void XTalk4sideHistogramManager::SetBinContentXTalk4side_horiz(double x, double y, double value_XTalk, double value_LY, double value_LY_XTalk){
    //---------2D----------//
    int x_binXly_2D = histogramsXTalk_2D_horiz_->GetXaxis()->FindBin(x);
    int y_binXly_2D = histogramsXTalk_2D_horiz_->GetYaxis()->FindBin(y);
    auto x_currentValue_norm_2D = histogramsXTalk_2D_horiz_->GetBinContent(x_binXly_2D, y_binXly_2D);

    int binX_2D = histogramsXTalk_2D_horiz_->GetXaxis()->FindBin(x);
    int binY_2D = histogramsXTalk_2D_horiz_->GetYaxis()->FindBin(y);
    auto currentCountX_2D = histogramsXTalk_2D_horiz_->GetBinContent(binX_2D, binY_2D);

    double x_newAvr_norm_2D = (x_currentValue_norm_2D * currentCountX_2D + value_XTalk)/(currentCountX_2D + 1);

    averageXTalk_horiz_->Fill(value_XTalk);
    dist_LY_XTalk_horiz_->Fill(value_LY_XTalk);
    dist_LY_horiz_->Fill(value_LY);
    pixelsData_horiz_LY_cube_[x_binXly_2D][y_binXly_2D].push_back(value_XTalk);
}

void XTalk4sideHistogramManager::SetBinContentXTalk4side_vert(double x, double y, double value_XTalk, double value_LY, double value_LY_XTalk){
    //---------2D----------//
    int x_binXly_2D = histogramsXTalk_2D_vert_->GetXaxis()->FindBin(x);
    int y_binXly_2D = histogramsXTalk_2D_vert_->GetYaxis()->FindBin(y);
    auto x_currentValue_norm_2D = histogramsXTalk_2D_vert_->GetBinContent(x_binXly_2D, y_binXly_2D);

    int binX_2D = histogramsXTalk_2D_vert_->GetXaxis()->FindBin(x);
    int binY_2D = histogramsXTalk_2D_vert_->GetYaxis()->FindBin(y);
    auto currentCountX_2D = histogramsXTalk_2D_vert_->GetBinContent(binX_2D, binY_2D);

    double x_newAvr_norm_2D = (x_currentValue_norm_2D * currentCountX_2D + value_XTalk)/(currentCountX_2D + 1);

    averageXTalk_vert_->Fill(value_XTalk);
    dist_LY_XTalk_vert_->Fill(value_LY_XTalk);
    dist_LY_vert_->Fill(value_LY);
    pixelsData_vert_LY_cube_[x_binXly_2D][y_binXly_2D].push_back(value_XTalk);
}

void XTalk4sideHistogramManager::WriteHistograms(TDirectory* dirName){
    auto average_XTalk_2D = dirName->mkdir("average_XTalk_2D");
    average_XTalk_2D->cd();
    histogramsXTalk_2D_up_->Draw();
    histogramsXTalk_2D_up_->Write();

    histogramsXTalk_2D_down_->Draw();
    histogramsXTalk_2D_down_->Write();

    histogramsXTalk_2D_right_->Draw();
    histogramsXTalk_2D_right_->Write();

    histogramsXTalk_2D_left_->Draw();
    histogramsXTalk_2D_left_->Write();

    histogramsXTalk_2D_horiz_->Draw();
    histogramsXTalk_2D_horiz_->Write();

    histogramsXTalk_2D_vert_->Draw();
    histogramsXTalk_2D_vert_->Write();

    auto average_XTalk_for_4_side = dirName->mkdir("average_XTalk_for_4_side");
    average_XTalk_for_4_side->cd();

//    averageXTalk_right_->Fit("landau", "r", " ", 0.02, 0.5);
    averageXTalk_right_->Draw();
    averageXTalk_right_->Write();

//    averageXTalk_left_->Fit("landau", "r", " ", 0.03, 0.5);
    averageXTalk_left_->Draw();
    averageXTalk_left_->Write();

//    averageXTalk_up_->Fit("landau", "r", " ", 0.02, 0.4);
    averageXTalk_up_->Draw();
    averageXTalk_up_->Write();

//    averageXTalk_down_->Fit("landau", "r", " ", 0.02, 0.4);
    averageXTalk_down_->Draw();
    averageXTalk_down_->Write();

    averageXTalk_horiz_->Draw();
    averageXTalk_horiz_->Write();

    averageXTalk_vert_->Draw();
    averageXTalk_vert_->Write();

    auto LY_XTalk_for_4_side = dirName->mkdir("LY_XTalk_for_4_side");
    LY_XTalk_for_4_side->cd();

    TF1* fit_gauss = new TF1("gaus", "gaus", 15, 40); // новые каты на световыход

    dist_LY_XTalk_up_->Draw();
    dist_LY_XTalk_up_->Write();

    dist_LY_XTalk_down_->Draw();
    dist_LY_XTalk_down_->Write();

    dist_LY_XTalk_right_->Draw();
    dist_LY_XTalk_right_->Write();

    dist_LY_XTalk_left_->Draw();
    dist_LY_XTalk_left_->Write();

    dist_LY_XTalk_horiz_->Draw();
    dist_LY_XTalk_horiz_->Write();

    dist_LY_XTalk_vert_->Draw();
    dist_LY_XTalk_vert_->Write();

    auto LY_for_4_side = dirName->mkdir("LY_for_4_side");
    LY_for_4_side->cd();

    dist_LY_up_->Fit( fit_gauss, "qr", " ", 15, 35);
    dist_LY_up_->Draw();
    dist_LY_up_->Write();

    dist_LY_down_->Fit( fit_gauss, "qr", " ", 15, 35);
    dist_LY_down_->Draw();
    dist_LY_down_->Write();

    dist_LY_right_->Fit( fit_gauss, "qr", " ", 15, 35);
    dist_LY_right_->Draw();
    dist_LY_right_->Write();

//    dist_LY_left_->Fit("landau", "r", " ", 1, 7);
    dist_LY_left_->Fit( fit_gauss, "qr", " ", 15, 35);
    dist_LY_left_->Draw();
    dist_LY_left_->Write();

    dist_LY_horiz_->Fit( fit_gauss, "qr", " ", 15, 35);
    dist_LY_horiz_->Draw();
    dist_LY_horiz_->Write();

    dist_LY_vert_->Fit( fit_gauss, "qr", " ", 15, 35);
    dist_LY_vert_->Draw();
    dist_LY_vert_->Write();

    auto all_2D = dirName->mkdir(Form("all_2D"));
    auto right_2D = all_2D->mkdir(Form("right_2D"));
    auto left_2D = all_2D->mkdir(Form("left_2D"));
    auto up_2D = all_2D->mkdir(Form("up_2D"));
    auto down_2D = all_2D->mkdir(Form("down_2D"));
    auto horiz_2D = all_2D->mkdir(Form("horiz_2D"));
    auto vert_2D = all_2D->mkdir(Form("vert_2D"));

    for (int x = 0; x < CUBES_SLICE_X; ++x) {
        for (int y = 0; y < CUBES_SLICE_Y; ++y) {
            if (!pixelsData_right_LY_cube_[x][y].empty()) {
                this->FillPixelData(right_2D, &pixelsData_right_LY_cube_[x][y], histogramsXTalk_2D_right_.get(),
                                    Form("pixel_%i_%i_right", x, y), x, y,0.5);
            }
            if (!pixelsData_left_LY_cube_[x][y].empty()) {
                this->FillPixelData(left_2D, &pixelsData_left_LY_cube_[x][y], histogramsXTalk_2D_left_.get(),
                                    Form("pixel_%i_%i_left", x, y), x, y,0.5);
            }
            if (!pixelsData_up_LY_cube_[x][y].empty()) {
                this->FillPixelData(up_2D, &pixelsData_up_LY_cube_[x][y], histogramsXTalk_2D_up_.get(),
                                    Form("pixel_%i_%i_up", x, y), x, y,0.5);
            }
            if (!pixelsData_down_LY_cube_[x][y].empty()) {
                this->FillPixelData(down_2D, &pixelsData_down_LY_cube_[x][y], histogramsXTalk_2D_down_.get(),
                                    Form("pixel_%i_%i_down", x, y), x, y, 0.5);
            }
            if (!pixelsData_horiz_LY_cube_[x][y].empty()) {
                this->FillPixelData(horiz_2D, &pixelsData_horiz_LY_cube_[x][y], histogramsXTalk_2D_horiz_.get(),
                                    Form("pixel_%i_%i_horiz", x, y), x, y, 0.5);
            }
            if (!pixelsData_vert_LY_cube_[x][y].empty()) {
                this->FillPixelData(horiz_2D, &pixelsData_vert_LY_cube_[x][y], histogramsXTalk_2D_vert_.get(),
                                    Form("pixel_%i_%i_vert", x, y), x, y, 0.5);
            }
        }
    }
    all_2D->cd();
    histogramsXTalk_2D_right_->Write();
    histogramsXTalk_2D_left_->Write();
    histogramsXTalk_2D_up_->Write();
    histogramsXTalk_2D_down_->Write();
    histogramsXTalk_2D_horiz_->Write();
    histogramsXTalk_2D_vert_->Write();

    auto all_1D = dirName->mkdir(Form("all_1D"));
    auto right_1D = all_1D->mkdir(Form("right_1D"));
    auto left_1D = all_1D->mkdir(Form("left_1D"));
    auto up_1D = all_1D->mkdir(Form("up_1D"));
    auto down_1D = all_1D->mkdir(Form("down_1D"));
    auto horiz_1D = all_1D->mkdir(Form("horiz_1D"));
    auto vert_1D = all_1D->mkdir(Form("vert_1D"));

    this->FillPixelAverageData_Left_Right(right_1D, pixelsData_right_LY_cube_, pixelsAverageXTalk_right_.get(),
                                          "pixel_average_right", 0.5);
    this->FillPixelAverageData_Left_Right(left_1D, pixelsData_left_LY_cube_, pixelsAverageXTalk_left_.get(),
                                          "pixel_average_left", 0.5);
    this->FillPixelAverageData_Up_Down(up_1D, pixelsData_up_LY_cube_, pixelsAverageXTalk_up_.get(),
                                       "pixel_average_up", 0.5);
    this->FillPixelAverageData_Up_Down(down_1D, pixelsData_down_LY_cube_, pixelsAverageXTalk_down_.get(),
                                       "pixel_average_down", 0.5);

    this->FillPixelAverageData_Horiz(horiz_1D, pixelsData_horiz_LY_cube_, pixelsAverageXTalk_horiz_.get(),
                                       "pixel_average_horiz", 0.5);
    this->FillPixelAverageData_Vert(vert_1D, pixelsData_vert_LY_cube_, pixelsAverageXTalk_vert_.get(),
                                     "pixel_average_vert", 0.5);

    all_1D->cd();
    pixelsAverageXTalk_right_->GetYaxis()->SetRangeUser(0, 0.1);
    pixelsAverageXTalk_right_->Write();
    pixelsAverageXTalk_left_->GetYaxis()->SetRangeUser(0, 0.1);
    pixelsAverageXTalk_left_->Write();
    pixelsAverageXTalk_up_->GetYaxis()->SetRangeUser(0, 0.1);
    pixelsAverageXTalk_up_->Write();
    pixelsAverageXTalk_down_->GetYaxis()->SetRangeUser(0, 0.1);
    pixelsAverageXTalk_down_->Write();
    pixelsAverageXTalk_horiz_->GetYaxis()->SetRangeUser(0, 0.1);
    pixelsAverageXTalk_horiz_->Write();
    pixelsAverageXTalk_vert_->GetYaxis()->SetRangeUser(0, 0.1);
    pixelsAverageXTalk_vert_->Write();
}

void XTalk4sideHistogramManager::FillPixelData(TDirectory* dir, std::vector<double>* data, TH2F* hist, const char* title, const int x, const int y, const double cut) {
    dir->cd();
    auto tempHist =  new TH1F(title, title, XTALK_BINNING);
    for (const double & it : *data) {
        tempHist->Fill(it);
    }

    if (tempHist->GetEntries() > 500) { // спец. кат, было 25000 для рана 125-128
        double pixel_mean = tempHist->GetMean();
        tempHist->Draw();
        tempHist->GetXaxis()->SetTitle("Crosstalk");
        tempHist->GetYaxis()->SetTitle("N");
        tempHist->Write();
        if ((pixel_mean == 0 || pixel_mean > 0) && pixel_mean < cut) {
            hist->SetBinContent(x, y, pixel_mean);
        }
    }
    tempHist->Delete();
}

void XTalk4sideHistogramManager::FillPixelAverageData_Left_Right(TDirectory* dir, const std::vector<double>(&data)[CUBES_SLICE_X][CUBES_SLICE_Y], TH1F* hist, const std::string& title, double cut) {
    dir->cd();
    std::unique_ptr<TH1F> tempHist[CUBES_SLICE_X];
    for (int x = 0; x < CUBES_SLICE_X; ++x) {
        std::string name = title + '_' + std::to_string(x);
        tempHist[x] =  std::make_unique<TH1F>(name.c_str(), name.c_str(), XTALK_BINNING);
    }
    for (int x = 0; x < CUBES_SLICE_X; ++x) {
        for (int y = 0; y < CUBES_SLICE_Y; ++y) {
            for (const auto it : data[x][y]) {
                tempHist[x]->Fill(it);
            }
        }

        if (tempHist[x]->GetEntries() > 500) { // спец. кат
            tempHist[x]->Draw();
            tempHist[x]->GetXaxis()->SetTitle("Crosstalk");
            tempHist[x]->GetYaxis()->SetTitle("N");
            tempHist[x]->Write();

            double pixel_mean = tempHist[x]->GetMean();
            if ((pixel_mean == 0 ||  pixel_mean > 0) && pixel_mean < cut) {
                hist->SetBinContent(x, pixel_mean);
            }
        }
    }
}

void XTalk4sideHistogramManager::FillPixelAverageData_Up_Down(TDirectory* dir, const std::vector<double>(&data)[CUBES_SLICE_X][CUBES_SLICE_Y], TH1F* hist, const std::string& title, double cut) {
    dir->cd();
    std::unique_ptr<TH1F> tempHist[CUBES_SLICE_Y];
    for (int y = 0; y < CUBES_SLICE_Y; ++y) {
        std::string name = title + '_' + std::to_string(y);
        tempHist[y] =  std::make_unique<TH1F>(name.c_str(), name.c_str(), XTALK_BINNING);
    }
    for (int y = 0; y < CUBES_SLICE_Y; ++y) {
        for (int x = 0; x < CUBES_SLICE_X; ++x) {
            for (const auto it : data[x][y]) {
                tempHist[y]->Fill(it);
            }
        }
        if ( tempHist[y]->GetEntries() > 500) { // спец. кат
            tempHist[y]->Draw();
            tempHist[y]->GetXaxis()->SetTitle("Crosstalk");
            tempHist[y]->GetYaxis()->SetTitle("N");
            tempHist[y]->Write();

            double pixel_mean = tempHist[y]->GetMean();
            if (pixel_mean > 0 && pixel_mean < cut) {
                hist->SetBinContent(y, pixel_mean);
            }
        }
    }
}

void XTalk4sideHistogramManager::FillPixelAverageData_Horiz(TDirectory* dir, const std::vector<double>(&data)[CUBES_SLICE_X][CUBES_SLICE_Y], TH1F* hist, const std::string& title, double cut) {
    dir->cd();
    std::unique_ptr<TH1F> tempHist[CUBES_SLICE_X];
    for (int x = 0; x < CUBES_SLICE_X; ++x) {
        std::string name = title + '_' + std::to_string(x);
        tempHist[x] =  std::make_unique<TH1F>(name.c_str(), name.c_str(), XTALK_BINNING);
    }
    for (int x = 0; x < CUBES_SLICE_X; ++x) {
        for (int y = 0; y < CUBES_SLICE_Y; ++y) {
            for (const auto it : data[x][y]) {
                tempHist[x]->Fill(it);
            }
        }

        if (tempHist[x]->GetEntries() > 500) { // спец. кат
            tempHist[x]->Draw();
            tempHist[x]->GetXaxis()->SetTitle("Crosstalk");
            tempHist[x]->GetYaxis()->SetTitle("N");
            tempHist[x]->Write();

            double pixel_mean = tempHist[x]->GetMean();
            if (pixel_mean > 0 && pixel_mean < cut) {
                hist->SetBinContent(x, pixel_mean);
            }
        }
    }
}

void XTalk4sideHistogramManager::FillPixelAverageData_Vert(TDirectory* dir, const std::vector<double>(&data)[CUBES_SLICE_X][CUBES_SLICE_Y], TH1F* hist, const std::string& title, double cut) {
        dir->cd();
        std::unique_ptr<TH1F> tempHist[CUBES_SLICE_Y];
        for (int y = 0; y < CUBES_SLICE_Y; ++y) {
            std::string name = title + '_' + std::to_string(y);
            tempHist[y] =  std::make_unique<TH1F>(name.c_str(), name.c_str(), XTALK_BINNING);
        }
        for (int y = 0; y < CUBES_SLICE_Y; ++y) {
            for (int x = 0; x < CUBES_SLICE_X; ++x) {
                for (const auto it : data[x][y]) {
                    tempHist[y]->Fill(it);
                }
            }
            if ( tempHist[y]->GetEntries() > 500) { // спец. кат, раньше было 25000
                tempHist[y]->Draw();
                tempHist[y]->GetXaxis()->SetTitle("Crosstalk");
                tempHist[y]->GetYaxis()->SetTitle("N");
                tempHist[y]->Write();

                double pixel_mean = tempHist[y]->GetMean();
                if (pixel_mean > 0 && pixel_mean < cut) {
                    hist->SetBinContent(y, pixel_mean);
                }
            }
        }
    }