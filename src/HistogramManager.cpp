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

        currentThreshold = maxContent * 0.35; // кат по количеству событий, можно поставить и 40

        if (cube == 92) {
            currentThreshold = maxContent * 0.45; // проблемный куб, нефизичный кат
        }

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

void LyHistogramManager::FillAveragePlots(const double x, const double y, const double x_ly, const double y_ly) {
    const int binX = histogramsAverageLY_gaus_layer_->GetXaxis()->FindBin(x);
    const int binY = histogramsAverageLY_gaus_layer_->GetYaxis()->FindBin(y);

    pixelsAverageSummLY_[binX][binY].push_back(x_ly+y_ly);
    pixelsAverageXLY_[binX][binY].push_back(x_ly);
    pixelsAverageYLY_[binX][binY].push_back(y_ly);
}

void LyHistogramManager::NormalizeCubes() const {
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

void LyHistogramManager::NormalizeFibers() const{
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
                        if (ch <= 20){ // z волокна идут по пучку
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
    // this->NormalizeCubes();
     this->NormalizeFibers();
//    averageCubeResponse_->Write();
//    averageCount_->Write();
}

void LyHistogramManager::FillPixelData(TDirectory* dir, std::vector<double>* data, TF1* fit, TH2F* hist, const char* title, const int x, const int y, const double cut) {
    dir->cd();
    const auto tempHist =  new TH1F(title, title,GAUS_BINNING);
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
        tempHist->Draw();
        tempHist->GetXaxis()->SetTitle("LY [p.e.]");
        tempHist->GetYaxis()->SetTitle("Number");
        tempHist->Write();

        double pixelLY = fit->GetParameter(1);
        double pixelSigma = fit->GetParameter(2);
        double pixelSigmaErr = fit->GetParError(2);

        double chi2 = fit->GetChisquare();
        int ndf = fit->GetNDF();

        if ((pixelLY > 0 && pixelLY < cut) && ndf != 0 && (pixelSigmaErr / pixelSigma < 0.5)){ // кат на отрезание плохого фитирования
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

void LyHistogramManager::PreparePlotsForGeometry(TDirectory* dirName, const std::map<int,CubePosition>* locationMap) {
    auto Geometry = dirName->mkdir("PlotsForGeometry");
    std::vector<CubeMapLY> x_avr;
    std::vector<CubeMapLY> y_avr;
    x_avr.reserve(DEFAULT_BINNING);
    y_avr.reserve(DEFAULT_BINNING);

    std::array<std::vector<CubeMapLY>,CUBES_COUNT> ly_x_on_x{};
    std::array<std::vector<CubeMapLY>,CUBES_COUNT> ly_y_on_y{};

    auto x_LY_on_X_projection = std::make_unique<TH1F>("x_LY_on_X_projection_AVR", "x_LY_on_X_projection_AVR", SINGLE_CUBE_BINNING_1D);
    std::array<std::unique_ptr<TH1F>,DEFAULT_BINNING> x_LY_on_X_bin;
    x_LY_on_X_projection.get()->GetXaxis()->SetTitle("X[mm]");
    x_LY_on_X_projection.get()->GetYaxis()->SetTitle("LY[p.e.]");

    auto y_LY_on_Y_projection = std::make_unique<TH1F>("y_LY_on_Y_projection_AVR", "y_LY_on_Y_projection_AVR", SINGLE_CUBE_BINNING_1D);
    std::array<std::unique_ptr<TH1F>,DEFAULT_BINNING> y_LY_on_Y_bin;
    y_LY_on_Y_projection.get()->GetXaxis()->SetTitle("Y[mm]");
    y_LY_on_Y_projection.get()->GetYaxis()->SetTitle("LY[p.e.]");

    for (int bin = 0; bin < DEFAULT_BINNING; ++bin) {
        x_LY_on_X_bin[bin] = std::make_unique<TH1F>(Form("x_LY_on_X_AVR_bin_%d",bin), Form("x_LY_on_X_AVR_bin_%d",bin), GAUS_BINNING);
        y_LY_on_Y_bin[bin] = std::make_unique<TH1F>(Form("y_LY_on_Y_AVR_bin_%d",bin), Form("y_LY_on_Y_AVR_bin_%d",bin), GAUS_BINNING);
        x_LY_on_X_bin[bin].get()->GetXaxis()->SetTitle("LY[p.e.]");
        x_LY_on_X_bin[bin].get()->GetYaxis()->SetTitle("Number");
        y_LY_on_Y_bin[bin].get()->GetXaxis()->SetTitle("LY[p.e.]");
        y_LY_on_Y_bin[bin].get()->GetYaxis()->SetTitle("Number");
    }
    std::array<std::unique_ptr<TH1F>, CUBES_COUNT> x_LY_X_projection;
    std::array<std::unique_ptr<TH1F>, CUBES_COUNT> y_LY_Y_projection;

    for (int cube = 0; cube <CUBES_COUNT; ++cube) {
        int layer = (cube/3)%3;

        ly_x_on_x[cube].reserve(DEFAULT_BINNING);
        ly_y_on_y[cube].reserve(DEFAULT_BINNING);

        if (locationMap->find(cube) == locationMap->end()) {
            std::cerr << "Warning: Cube " << cube << " not found in locationMap. Skipping..." << std::endl;
            continue;
        }

        const auto cube_position = locationMap->at(cube);
        auto cube_plots = Geometry->mkdir(Form("cube_%d", cube));
        cube_plots->cd();
        {
            auto binX_min = histogramsLY_normalized_layer_[0]->GetXaxis()->FindBin(cube_position._x_min);
            auto binX_max = histogramsLY_normalized_layer_[0]->GetXaxis()->FindBin(cube_position._x_max);
            auto binY_min = histogramsLY_normalized_layer_[0]->GetYaxis()->FindBin(cube_position._y_min);
            auto binY_max = histogramsLY_normalized_layer_[0]->GetYaxis()->FindBin(cube_position._y_max);

            auto summ_LY_X_projection = std::make_unique<TH1F>(Form("summ_LY_X_projection_%d", cube), Form("summ_LY_X_projection_%d", cube), SINGLE_CUBE_BINNING_1D);
            auto summ_LY_Y_projection = std::make_unique<TH1F>(Form("summ_LY_Y_projection_%d", cube), Form("summ_LY_Y_projection_%d", cube), SINGLE_CUBE_BINNING_1D);
            summ_LY_X_projection.get()->GetXaxis()->SetTitle("X[mm]");
            summ_LY_X_projection.get()->GetYaxis()->SetTitle("LY[p.e.]");
            summ_LY_Y_projection.get()->GetXaxis()->SetTitle("Y[mm]");
            summ_LY_Y_projection.get()->GetYaxis()->SetTitle("LY[p.e.]");

            x_LY_X_projection[cube] = std::make_unique<TH1F>(Form("x_LY_X_projection_%d", cube), Form("x_LY_X_projection_%d", cube), SINGLE_CUBE_BINNING_1D);
            auto x_LY_Y_projection = std::make_unique<TH1F>(Form("x_LY_Y_projection_%d", cube), Form("x_LY_Y_projection_%d", cube), SINGLE_CUBE_BINNING_1D);
            x_LY_X_projection[cube].get()->GetXaxis()->SetTitle("X[mm]");
            x_LY_X_projection[cube].get()->GetYaxis()->SetTitle("LY[p.e.]");
            x_LY_Y_projection.get()->GetXaxis()->SetTitle("Y[mm]");
            x_LY_Y_projection.get()->GetYaxis()->SetTitle("LY[p.e.]");

            auto y_LY_X_projection = std::make_unique<TH1F>(Form("y_LY_X_projection_%d", cube), Form("y_LY_X_projection_%d", cube), SINGLE_CUBE_BINNING_1D);
            y_LY_Y_projection[cube] = std::make_unique<TH1F>(Form("y_LY_Y_projection_%d", cube), Form("y_LY_Y_projection_%d", cube), SINGLE_CUBE_BINNING_1D);
            y_LY_X_projection.get()->GetXaxis()->SetTitle("X[mm]");
            y_LY_X_projection.get()->GetYaxis()->SetTitle("LY[p.e.]");
            y_LY_Y_projection[cube].get()->GetXaxis()->SetTitle("Y[mm]");
            y_LY_Y_projection[cube].get()->GetYaxis()->SetTitle("LY[p.e.]");

            TF1* fit = new TF1("gaus","gaus", 10 , 200);
            for (int x = binX_min; x <= binX_max; ++x) {
                FillPixelProjectionXData(cube_plots, pixelsDataSummLY_, fit, summ_LY_X_projection.get(), Form("projection_layer_%i_for_x_%i_summ", layer, x ), x, binX_min, binX_max, 0, binY_min, binY_max, layer, 140);
                const auto ly_x = FillPixelProjectionXData(cube_plots, pixelsDataXLY_, fit, x_LY_X_projection[cube].get(), Form("projection_layer_%i_for_x_%i_x", layer, x ), x, binX_min, binX_max, 0, binY_min, binY_max, layer, 70);
                FillPixelProjectionXData(cube_plots, pixelsDataYLY_, fit, y_LY_X_projection.get(), Form("projection_layer_%i_for_x_%i_y", layer, x ), x, binX_min, binX_max, 0, binY_min, binY_max, layer, 70);
                ly_x_on_x[cube].push_back({x - binX_min,0,ly_x});
                for (int y = binY_min; y < binY_max && x - binX_min < DEFAULT_BINNING; ++y) {
                    for (const auto & ly: pixelsDataXLY_[layer][x][y]) {
                        x_LY_on_X_bin[x - binX_min].get()->Fill(ly);
                    }
                }
            }
            cube_plots->cd();
            summ_LY_X_projection->Write();
            x_LY_X_projection[cube]->Write();
            y_LY_X_projection->Write();
            for (int y = binY_min; y < binY_max; ++y) {
                FillPixelProjectionYData(cube_plots, pixelsDataSummLY_, fit, summ_LY_Y_projection.get(), Form("projection_layer_%i_for_y_%i_summ", layer, y ), 0, binX_min, binX_max, y, binY_min, binY_max, layer, 140);
                FillPixelProjectionYData(cube_plots, pixelsDataXLY_, fit, x_LY_Y_projection.get(), Form("projection_layer_%i_for_y_%i_x", layer, y ), 0, binX_min, binX_max, y, binY_min, binY_max, layer, 70);
                const auto ly_y = FillPixelProjectionYData(cube_plots, pixelsDataYLY_, fit, y_LY_Y_projection[cube].get(), Form("projection_layer_%i_for_y_%i_y", layer, y ), 0, binX_min, binX_max, y, binY_min, binY_max, layer, 70);
                ly_y_on_y[cube].push_back({0,y - binY_min,ly_y});
                for (int x = binX_min; x <= binX_max; ++x) {
                    for (const auto & ly: pixelsDataYLY_[layer][x][y]) {
                        y_LY_on_Y_bin[y - binY_min].get()->Fill(ly);
                    }
                }
            }
            cube_plots->cd();
            summ_LY_Y_projection->Write();
            x_LY_Y_projection->Write();
            y_LY_Y_projection[cube]->Write();
        }
    }

    auto avr_Cube = Geometry->mkdir("Avr_Cube");
    auto X = avr_Cube->mkdir("X");
    auto Y = avr_Cube->mkdir("Y");
    TF1* fit = new TF1("gaus","gaus", 10 , 100);
    for (int bin = 0; bin < DEFAULT_BINNING; ++bin) {
        std::unique_ptr<TSpectrum> tSpectrum = std::make_unique<TSpectrum>(15);
        int nFound = tSpectrum->Search(x_LY_on_X_bin[bin].get(),10,"nobackground new",0.1);
        {
            if(nFound > 0) {
                double *xpeaks = tSpectrum->GetPositionX();
                for (uint peak = 0; peak < nFound; ++peak) {
                    if (xpeaks[peak] > 15) {
                        x_LY_on_X_bin[bin]->Fit(fit,"qr+", "",xpeaks[peak] - 15, xpeaks[peak] + 15);
                        break;
                    }
                }
            }
            X->cd();
            x_LY_on_X_bin[bin]->Write();

            double pixelLY = fit->GetParameter(1);
            x_avr.push_back({bin,0,pixelLY});

            if (pixelLY > 0) {
                x_LY_on_X_projection->SetBinContent(bin, pixelLY);
            }
        }

        {
            nFound = tSpectrum->Search(y_LY_on_Y_bin[bin].get(),10,"nobackground new",0.1);
            if(nFound > 0) {
                double *xpeaks = tSpectrum->GetPositionX();
                for (uint peak = 0; peak < nFound; ++peak) {
                    if (xpeaks[peak] > 15) { // новый кат?
                        y_LY_on_Y_bin[bin]->Fit(fit,"qr+", "",xpeaks[peak] - 15, xpeaks[peak] + 15); // новый кат?
                        break;
                    }
                }
            }
            Y->cd();
            y_LY_on_Y_bin[bin]->Write();

            double pixelLY = fit->GetParameter(1);
            y_avr.push_back({0,bin,pixelLY});
            if (pixelLY > 0 ) {
                y_LY_on_Y_projection->SetBinContent(bin, pixelLY);
            }
        }
    }

    for (int cube = 0; cube < CUBES_COUNT; ++cube) {
        auto cube_plots = Geometry->GetDirectory(Form("cube_%d", cube));
        auto Diff_X_on_X_LY = std::make_unique<TH1F>(Form("Diff_X_on_X_LY_%d", cube), Form("Diff_X_on_X_LY_%d", cube), SINGLE_CUBE_BINNING_1D);
        auto Diff_Y_on_Y_LY = std::make_unique<TH1F>(Form("Diff_Y_on_Y_LY_%d", cube), Form("Diff_Y_on_Y_LY_%d", cube), SINGLE_CUBE_BINNING_1D);
        Diff_X_on_X_LY.get()->GetXaxis()->SetTitle("X[mm]");
        Diff_X_on_X_LY.get()->GetYaxis()->SetTitle("LY[p.e.]");
        Diff_Y_on_Y_LY.get()->GetXaxis()->SetTitle("Y[mm]");
        Diff_Y_on_Y_LY.get()->GetYaxis()->SetTitle("LY[p.e.]");

        TCanvas c1 = TCanvas(Form("Same_Y_on_Y_cube_%d_LY_avrLY",cube), Form("Same_Y_on_Y_cube_%d_LY_avrLY",cube), 800, 600);
        y_LY_on_Y_projection->SetTitle(Form("Same_Y_on_Y_cube_%d_LY_avrLY",cube));
        y_LY_on_Y_projection->Draw("HIST");
        y_LY_on_Y_projection->SetLineColor(kRed);
        y_LY_Y_projection[cube]->Draw("HIST SAME");
        y_LY_Y_projection[cube]->SetLineColor(kBlue);
        c1.Update();

        TCanvas c2 = TCanvas(Form("Same_X_on_X_cube_%d_LY_avrLY",cube), Form("Same_X_on_X_cube_%d_LY_avrLY",cube), 800, 600);
        x_LY_on_X_projection->SetTitle(Form("Same_X_on_X_cube_%d_LY_avrLY",cube));
        x_LY_on_X_projection->Draw("HIST_2");
        x_LY_on_X_projection->SetLineColor(kRed);
        x_LY_X_projection[cube]->Draw("HIST_2 SAME");
        x_LY_X_projection[cube]->SetLineColor(kBlue);
        c2.Update();

        for (int bin = 0; bin < DEFAULT_BINNING; ++bin) {
            {
                auto it = std::find_if(ly_x_on_x[cube].begin(), ly_x_on_x[cube].end(), [&](const auto ly) {
                   return bin == ly.x_;
               });
                auto it_avr = std::find_if(x_avr.begin(), x_avr.end(), [&](const auto ly) {
                    return bin == ly.x_;
                });
                if (it != ly_x_on_x[cube].end() && it_avr != x_avr.end())
                    Diff_X_on_X_LY->SetBinContent(bin, (it->ly_ - it_avr->ly_)/it_avr->ly_);
            }
            {
                auto it = std::find_if(ly_y_on_y[cube].begin(), ly_y_on_y[cube].end(), [&](const auto ly) {
                   return bin == ly.y_;
               });
                auto it_avr = std::find_if(y_avr.begin(), y_avr.end(), [&](const auto ly) {
                    return bin == ly.y_;
                });
                if (it != ly_y_on_y[cube].end() && it_avr != y_avr.end())
                    Diff_Y_on_Y_LY->SetBinContent(bin, (it->ly_ - it_avr->ly_)/it_avr->ly_);
            }
        }
        cube_plots->cd();
        Diff_X_on_X_LY->Write();
        Diff_Y_on_Y_LY->Write();
        c1.Write();
        c2.Write();
    }

    avr_Cube->cd();
    x_LY_on_X_projection->Write();
    y_LY_on_Y_projection->Write();
}
void LyHistogramManager::PreparePlotsForMC(TDirectory* dirName, const std::map<int,CubePosition>* locationMap) {
    auto MC = dirName->mkdir("PlotsForMC");
    std::vector<CubeMapLY> ly_summ_avr;
    std::vector<CubeMapLY> ly_x_avr;
    std::vector<CubeMapLY> ly_y_avr;
    ly_summ_avr.reserve(DEFAULT_BINNING*DEFAULT_BINNING);
    ly_x_avr.reserve(DEFAULT_BINNING*DEFAULT_BINNING);
    ly_y_avr.reserve(DEFAULT_BINNING*DEFAULT_BINNING);

    std::array<std::vector<CubeMapLY>,CUBES_COUNT> ly_summ{};
    std::array<std::vector<CubeMapLY>,CUBES_COUNT> ly_x{};
    std::array<std::vector<CubeMapLY>,CUBES_COUNT> ly_y{};

    for (int cube = 0; cube <CUBES_COUNT; ++cube) {
        ly_summ[cube].reserve(DEFAULT_BINNING*DEFAULT_BINNING);
        ly_x[cube].reserve(DEFAULT_BINNING*DEFAULT_BINNING);
        ly_y[cube].reserve(DEFAULT_BINNING*DEFAULT_BINNING);
        int layer = (cube/3)%3;

        if (locationMap->find(cube) == locationMap->end()) {
            std::cerr << "Warning: Cube " << cube << " not found in locationMap. Skipping..." << std::endl;
            continue;
        }
        const auto cube_position = locationMap->at(cube);


        if (histogramsLY_normalized_layer_.empty() || !histogramsLY_normalized_layer_[0]) {
            std::cerr << "Critical Error: histogramsLY_normalized_layer_ is empty!" << std::endl;
            return;
        }

        auto cube_plots = MC->mkdir(Form("cube_%d", cube));
        cube_plots->cd();
        auto cube_LY_sum_mean = std::make_unique<TH2F>(Form("cube_LY_Summ_Mean_%d", cube), Form("cube_LY_Summ_Mean_%d", cube), CUBES_BINNING);
        auto cube_LY_X_mean = std::make_unique<TH2F>(Form("cube_LY_X_Mean_%d", cube), Form("cube_LY_X_Mean_%d", cube), CUBES_BINNING);
        auto cube_LY_Y_mean = std::make_unique<TH2F>(Form("cube_LY_Y_Mean_%d", cube), Form("cube_LY_Y_Mean_%d", cube), CUBES_BINNING);
        cube_LY_sum_mean.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_sum_mean.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_sum_mean.get()->GetZaxis()->SetTitle("LY[p.e.]");
        cube_LY_X_mean.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_X_mean.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_X_mean.get()->GetZaxis()->SetTitle("LY[p.e.]");
        cube_LY_Y_mean.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_Y_mean.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_Y_mean.get()->GetZaxis()->SetTitle("LY[p.e.]");

        auto cube_LY_sum_error = std::make_unique<TH2F>(Form("cube_LY_Summ_Error_%d", cube), Form("cube_LY_Summ_Error_%d", cube), CUBES_BINNING);
        auto cube_LY_X_error = std::make_unique<TH2F>(Form("cube_LY_X_Error_%d", cube), Form("cube_LY_X_Error_%d", cube), CUBES_BINNING);
        auto cube_LY_Y_error = std::make_unique<TH2F>(Form("cube_LY_Y_Error_%d", cube), Form("cube_LY_Y_Error_%d", cube), CUBES_BINNING);
        cube_LY_sum_error.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_sum_error.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_sum_error.get()->GetZaxis()->SetTitle("Error[p.e.]");
        cube_LY_X_error.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_X_error.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_X_error.get()->GetZaxis()->SetTitle("Error[p.e.]");
        cube_LY_Y_error.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_Y_error.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_Y_error.get()->GetZaxis()->SetTitle("Error[p.e.]");

        auto cube_LY_sum_sigma = std::make_unique<TH2F>(Form("cube_LY_Summ_Sigma_%d", cube), Form("cube_LY_Summ_Sigma_%d", cube), CUBES_BINNING);
        auto cube_LY_X_sigma = std::make_unique<TH2F>(Form("cube_LY_X_Sigma_%d", cube), Form("cube_LY_X_Sigma_%d", cube), CUBES_BINNING);
        auto cube_LY_Y_sigma = std::make_unique<TH2F>(Form("cube_LY_Y_Sigma_%d", cube), Form("cube_LY_Y_Sigma_%d", cube), CUBES_BINNING);
        cube_LY_sum_sigma.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_sum_sigma.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_sum_sigma.get()->GetZaxis()->SetTitle("Sigma[p.e.]");
        cube_LY_X_sigma.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_X_sigma.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_X_sigma.get()->GetZaxis()->SetTitle("Sigma[p.e.]");
        cube_LY_Y_sigma.get()->GetXaxis()->SetTitle("X[mm]");
        cube_LY_Y_sigma.get()->GetYaxis()->SetTitle("Y[mm]");
        cube_LY_Y_sigma.get()->GetZaxis()->SetTitle("Sigma[p.e.]");

        {
            auto binX_min = histogramsLY_normalized_layer_[0]->GetXaxis()->FindBin(cube_position._x_min);
            auto binX_max = histogramsLY_normalized_layer_[0]->GetXaxis()->FindBin(cube_position._x_max);
            auto binY_min = histogramsLY_normalized_layer_[0]->GetYaxis()->FindBin(cube_position._y_min);
            auto binY_max = histogramsLY_normalized_layer_[0]->GetYaxis()->FindBin(cube_position._y_max);
            // std::cout << cube << " : " << AS_KV(binX_min) << ' ' << AS_KV(binX_max) << ' ' << AS_KV(binY_min) << ' ' << AS_KV(binY_max) <<std::endl;
            for (int x = binX_min; x <= binX_max; ++x) {
                for (int y = binY_min; y < binY_max; ++y) {
                    TF1* fit = new TF1("gaus","gaus", 10 , 200);
                    if (!pixelsDataSummLY_[layer][x][y].empty()) {
                        auto ly = FillPixelData(&pixelsDataSummLY_[layer][x][y], fit, cube_LY_sum_mean.get(), cube_LY_sum_error.get(), cube_LY_sum_sigma.get(),
                            Form("pixel_%i_%i_%i_summ", layer, x, y), x, y, 140);
                        if (ly > 0) {
                            ly_summ[cube].push_back({x - binX_min, y - binY_min, ly});
                        }
                    }
                    //----------------------------------------//----------------------------------------//
                    if (!pixelsDataXLY_[layer][x][y].empty()) {
                        fit = new TF1("gaus","gaus", 20 , 80);
                        auto ly = FillPixelData(&pixelsDataXLY_[layer][x][y], fit, cube_LY_X_mean.get(), cube_LY_X_error.get(), cube_LY_X_sigma.get(),
                            Form("pixel_%i_%i_%i_x",layer,x,y),x,y, 70);
                        if (ly > 0) {
                            ly_x[cube].push_back({x - binX_min, y - binY_min, ly});
                        }
                    }
                    //----------------------------------------//----------------------------------------//
                    if (!pixelsDataYLY_[layer][x][y].empty()) {
                        fit = new TF1("gaus","gaus", 20 , 80);
                        auto ly = FillPixelData(&pixelsDataYLY_[layer][x][y], fit, cube_LY_Y_mean.get(),cube_LY_Y_error.get(),cube_LY_Y_sigma.get(),
                       Form("pixel_%i_%i_%i_y",layer,x,y),x,y,70);
                        if (ly > 0) {
                            ly_y[cube].push_back({x - binX_min, y - binY_min, ly});
                        }
                    }
                    fit->Delete();
                }
            }
        }
        cube_LY_sum_mean->Write();
        cube_LY_X_mean->Write();
        cube_LY_Y_mean->Write();
        cube_LY_sum_error->Write();
        cube_LY_X_error->Write();
        cube_LY_Y_error->Write();
        cube_LY_sum_sigma->Write();
        cube_LY_X_sigma->Write();
        cube_LY_Y_sigma->Write();
    }
    auto cube_LY_sum_arv = std::make_unique<TH2F>("cube_LY_Summ_Mean_Avr", "cube_LY_Summ_Mean_Avr", SINGLE_CUBE_BINNING);
    auto cube_LY_X_arv = std::make_unique<TH2F>("cube_LY_X_Mean_Avr", "cube_LY_X_Mean_Avr", SINGLE_CUBE_BINNING);
    auto cube_LY_Y_arv = std::make_unique<TH2F>("cube_LY_Y_Mean_Avr", "cube_LY_Y_Mean_Avr", SINGLE_CUBE_BINNING);
    cube_LY_sum_arv->GetXaxis()->SetTitle("X[mm]");
    cube_LY_sum_arv->GetYaxis()->SetTitle("Y[mm]");
    cube_LY_sum_arv.get()->GetZaxis()->SetTitle("LY_X+Y[p.e.]");
    cube_LY_X_arv.get()->GetXaxis()->SetTitle("X[mm]");
    cube_LY_X_arv.get()->GetYaxis()->SetTitle("Y[mm]");
    cube_LY_X_arv.get()->GetZaxis()->SetTitle("LY_X[p.e.]");
    cube_LY_Y_arv.get()->GetXaxis()->SetTitle("X[mm]");
    cube_LY_Y_arv.get()->GetYaxis()->SetTitle("Y[mm]");
    cube_LY_Y_arv.get()->GetZaxis()->SetTitle("LY_Y[p.e.]");
    for (int x = 0; x < DEFAULT_BINNING; ++x) {
        for (int y = 0; y < DEFAULT_BINNING; ++y) {
            size_t count_summ = 0, count_x = 0, count_y = 0;
            CubeMapLY temp_summ{x,y,0};
            CubeMapLY temp_x{x,y,0};
            CubeMapLY temp_y{x,y,0};
            for (int cube = 0; cube < CUBES_COUNT; ++cube) {
                {
                    auto it = std::find_if(ly_summ[cube].begin(), ly_summ[cube].end(), [&](const auto ly) {
                       return x == ly.x_ && y == ly.y_;
                   });
                    if (it != ly_summ[cube].end()) {
                        temp_summ.ly_ += it->ly_;
                        ++count_summ;
                    }
                }
                {
                    auto it = std::find_if(ly_x[cube].begin(), ly_x[cube].end(), [&](const auto ly) {
                        return x == ly.x_ && y == ly.y_;
                    });
                    if (it != ly_x[cube].end()) {
                        temp_x.ly_ += it->ly_;
                        ++count_x;
                    }
                }
                {
                    auto it = std::find_if(ly_y[cube].begin(), ly_y[cube].end(), [&](const auto ly) {
                        return x == ly.x_ && y == ly.y_;
                    });
                    if (it != ly_y[cube].end()) {
                        temp_y.ly_ += it->ly_;
                        ++count_y;
                    }
                }
            }
            temp_summ.ly_ = temp_summ.ly_/count_summ;
            if (temp_summ.ly_ > 0) {
                cube_LY_sum_arv->SetBinContent(x,y,temp_summ.ly_);
                ly_summ_avr.push_back(std::move(temp_summ));
            }

            temp_x.ly_ = temp_x.ly_ / count_x;
            if (temp_x.ly_ > 0) {
                cube_LY_X_arv->SetBinContent(x,y,temp_x.ly_);
                ly_x_avr.push_back(std::move(temp_x));
            }

            temp_y.ly_ = temp_y.ly_ / count_y;
            if (temp_y.ly_ > 0) {
                cube_LY_Y_arv->SetBinContent(x,y,temp_y.ly_);
                ly_y_avr.push_back(std::move(temp_y));
            }
        }
    }
    MC->cd();
    cube_LY_sum_arv->Write();
    cube_LY_X_arv->Write();
    cube_LY_Y_arv->Write();

    for (int cube = 0; cube < CUBES_COUNT; ++cube) {
        auto cube_plots = MC->GetDirectory(Form("cube_%d", cube));
        cube_plots->cd();
        auto cube_LY_sum_diff = std::make_unique<TH2F>(Form("cube_LY_Summ_Diff_%d", cube), Form("cube_LY_Summ_Diff_%d", cube), SINGLE_CUBE_BINNING);
        auto cube_LY_X_diff = std::make_unique<TH2F>(Form("cube_LY_X_Diff_%d", cube), Form("cube_LY_X_Diff_%d", cube), SINGLE_CUBE_BINNING);
        auto cube_LY_Y_diff = std::make_unique<TH2F>(Form("cube_LY_Y_Diff_%d", cube), Form("cube_LY_Y_Diff_%d", cube), SINGLE_CUBE_BINNING);
        for (int x = 0; x < DEFAULT_BINNING - 1; ++x) {
            for (int y = 0; y < DEFAULT_BINNING - 1; ++y) {
                {
                    auto it_avr = std::find_if(ly_summ_avr.begin(), ly_summ_avr.end(), [&](const auto ly) {
                       return x == ly.x_ && y == ly.y_;
                   });

                    auto it = std::find_if(ly_summ[cube].begin(), ly_summ[cube].end(), [&](const auto ly) {
                       return x == ly.x_ && y == ly.y_;
                   });

                    if (it != ly_summ[cube].end()) {
                        auto ly_diff = (it->ly_ - it_avr->ly_)/it_avr->ly_;
                        cube_LY_sum_diff->SetBinContent(x,y,ly_diff);
                    }
                }
                {
                    auto it_avr = std::find_if(ly_x_avr.begin(), ly_x_avr.end(), [&](const auto ly) {
                        return x == ly.x_ && y == ly.y_;
                    });

                    auto it = std::find_if(ly_x[cube].begin(), ly_x[cube].end(), [&](const auto ly) {
                        return x == ly.x_ && y == ly.y_;
                    });

                    if (it != ly_x[cube].end()) {
                        auto ly_diff = (it->ly_ - it_avr->ly_)/it_avr->ly_;
                        cube_LY_X_diff->SetBinContent(x,y,ly_diff);
                    }
                }
                {
                    auto it_avr = std::find_if(ly_y_avr.begin(), ly_y_avr.end(), [&](const auto ly) {
                        return x == ly.x_ && y == ly.y_;
                    });

                    auto it = std::find_if(ly_y[cube].begin(), ly_y[cube].end(), [&](const auto ly) {
                        return x == ly.x_ && y == ly.y_;
                    });

                    if (it != ly_y[cube].end()) {
                        auto ly_diff = (it->ly_ - it_avr->ly_)/it_avr->ly_;
                        cube_LY_Y_diff->SetBinContent(x,y,ly_diff);
                    }
                }
            }
        }
        cube_LY_sum_diff->Write();
        cube_LY_X_diff->Write();
        cube_LY_Y_diff->Write();
    }
}
void LyHistogramManager::FillPixelFiberData(TDirectory* dir, std::vector<double>* data, TH1F* hist, TH1F* histAvr, TH1F* histCompr, const char* title, int position) {
    dir->cd();
    TH1F* tempHist =  new TH1F(title, title,GAUS_FIBER_BINNING);
    unsigned int peak1_count = 0;
    unsigned int peak2_count = 0;

    for (const double & it : *data) {
        tempHist->Fill(it);
        it < 18 ? ++peak1_count : ++peak2_count;
    }

    std::unique_ptr<TSpectrum> tSpectrum = std::make_unique<TSpectrum>(15);
    int nFound = tSpectrum->Search(tempHist,10,"nobackground new",0.1);
    TF1* fit = new TF1("gaus","gaus");
    if(nFound > 0) {
        double *xpeaks = tSpectrum->GetPositionX();
        for (uint peak = 0; peak < nFound; ++peak) {
            if (xpeaks[peak] >= 18 && xpeaks[peak] < 70) {
                fit = new TF1("gaus","gaus");
                tempHist->Fit(fit,"qr+", "",xpeaks[peak] - 15, xpeaks[peak] + 15);
                break;
            }
            if (xpeaks[peak] > 0 && xpeaks[peak] < 18) {
                fit = new TF1("landau","landau");
                tempHist->Fit(fit,"qr+", "", xpeaks[peak] - 5, xpeaks[peak] + 5);
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
        histCompr->SetBinContent(position, static_cast<double>(peak1_count)/static_cast<double>(peak2_count));

    double pixelLY = fit->GetParameter(1);
    if (pixelLY > 0 && pixelLY < 80) {
        hist->SetBinContent(position,pixelLY);
    }

    histAvr->SetBinContent(position, tempHist->GetMean());
    tempHist->Delete();
}
