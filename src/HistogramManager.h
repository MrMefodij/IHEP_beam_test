//
// Created by amefodev on 03.11.2024.
//

#ifndef CUBESLY_HISTOGRAMMANAGER_H
#define CUBESLY_HISTOGRAMMANAGER_H
#include <connection_map.h>
#include <vector>
#include <memory>
#include <TH1F.h>
#include <TH2F.h>
#include <string>
#include <TFile.h>
#include <TDirectory.h>
#include <ostream>
#include "defines.h"

class HistogramManager {
public:
    virtual ~HistogramManager() = default;
    explicit HistogramManager(unsigned int count = CUBES_COUNT);
    TH2F* getHistogram(size_t index);
    void fillHistograms(unsigned int  index, double x, double y);
    void SetBinContent(unsigned int  index, double x, double y, double value);
    void SetTitlesHist(const std::string& titleX, const std::string& titleY, const std::string& titleZ);
    virtual void WriteHistograms(TDirectory* dirName);
    void FindHighDensityArea(int threshold, std::ofstream& fileOut);
protected:
    std::vector<std::unique_ptr<TH2F>> histograms_;
};

class LyHistogramManager : public HistogramManager {
public:
    explicit LyHistogramManager(unsigned int count = CUBES_COUNT);

    virtual void SetTitlesHistLY(const std::string& titleX, const std::string& titleY, const std::string& titleZ);
    void SetBinContentLY(unsigned int index, double x, double y, double x_value, double y_value);
    virtual void SetBinContentNormalizedLY(unsigned int index, double x, double y, double x_value, double y_value, double z_value);
    void SetBinContentNormalizedLY(unsigned int index, double x, double y, double x_value, double y_value);
    void SetBinContentNormalizedLY_layer(unsigned int index, double x, double y, double x_value, double y_value);
    void SetBinContentNormalizedLY_layer(unsigned int index, double x, double y, double x_value, double y_value, double z_value);
    void SetFiberLY_layer(double x_fiber, double y_fiber, double x_position, double y_position,  double x_ly, double y_ly);
    void FillEventTime(double value);
    void SetBinContentAvrCubeResponse(double x, double y, double value);
    void WriteHistograms(TDirectory* dirName) override;
    virtual void SetTitlesHistLY_Fiber_I(const std::string& titleX, const std::string& titleY) ;
    virtual void SetTitlesHistLY_Fiber_A(const std::string& titleX, const std::string& titleY);
    virtual void SetTitlesHistLY_Fiber_A_reconstructed(const std::string& titleX, const std::string& titleY);
    void SetFiberLY_I(unsigned int index, double x, double y);
    void SetFiberLY_A(unsigned int index, double x, double y);
    void SetFiberLY_A_recon(unsigned int index, double x, double y);
    void FillAveragePlots(double x, double y, double x_ly, double y_ly);
    virtual void FillFiberLY(unsigned int ch_x, unsigned int ch_y, double ly_x, double ly_y);
    void PreparePlotsForMC(TDirectory* dirName, const std::map<unsigned int,CubePosition>* locationMap);
    void PreparePlotsForGeometry(TDirectory* dirName, const std::map<unsigned int,CubePosition>* locationMap);
protected:
    std::vector<std::unique_ptr<TH2F>> histogramsLY_;
    std::vector<std::unique_ptr<TH2F>> x_HistogramsLY_;
    std::vector<std::unique_ptr<TH2F>> y_HistogramsLY_;

    std::vector<std::unique_ptr<TH2F>> histogramsLY_normalized_;
    std::vector<std::unique_ptr<TH2F>> x_HistogramsLY_normalized_;
    std::vector<std::unique_ptr<TH2F>> y_HistogramsLY_normalized_;

    std::vector<std::unique_ptr<TH1F>> x_CubeLY_I;
    std::vector<std::unique_ptr<TH1F>> y_CubeLY_I;
    std::vector<std::unique_ptr<TH1F>> x_CubeLY_A;
    std::vector<std::unique_ptr<TH1F>> y_CubeLY_A;
    std::vector<std::unique_ptr<TH1F>> x_CubeLY_A_recon;
    std::vector<std::unique_ptr<TH1F>> y_CubeLY_A_recon;

    std::unique_ptr<TH1F> averageLY_;
    std::unique_ptr<TH1F> averageLY_X_;
    std::unique_ptr<TH1F> averageLY_Y_;
    std::unique_ptr<TH1F> averageTime_;

    std::unique_ptr<TH2F> averageCubeResponse_;
    std::unique_ptr<TH2F> averageCount_;

    std::vector<std::unique_ptr<TH2F>> histograms_layer_;
    std::vector<std::unique_ptr<TH2F>> histogramsLY_normalized_layer_;
    std::vector<std::unique_ptr<TH2F>> x_HistogramsLY_normalized_layer_;
    std::vector<std::unique_ptr<TH2F>> y_HistogramsLY_normalized_layer_;

    std::vector<std::unique_ptr<TH1D>> histogramsLY_normalized_layer_slice_X_;
    std::vector<std::unique_ptr<TH1D>> histogramsLY_normalized_layer_slice_Y_;
    std::vector<std::unique_ptr<TH1D>> x_HistogramsLY_normalized_layer_slice_;
    std::vector<std::unique_ptr<TH1D>> y_HistogramsLY_normalized_layer_slice_;

    std::vector<std::unique_ptr<TH1F>> histograms_fiberLY_normalized_layer_;

    std::vector<double> x_fiber_pixel_LY_[CHANNELS_NUMBER][CUBES_SLICE_Y];
    std::vector<double> y_fiber_pixel_LY_[CHANNELS_NUMBER][CUBES_SLICE_X];
    std::vector<std::unique_ptr<TH1F>> x_fiberLY_normalized_position_;
    std::vector<std::unique_ptr<TH1F>> y_fiberLY_normalized_position_;
    std::vector<std::unique_ptr<TH1F>> x_Avr_Mean_fiberLY_normalized_position_;
    std::vector<std::unique_ptr<TH1F>> y_Avr_Mean_fiberLY_normalized_position_;

    std::vector<std::unique_ptr<TH1F>> x_fiber_peak_comparison_;
    std::vector<std::unique_ptr<TH1F>> y_fiber_peak_comparison_;

    std::vector<std::unique_ptr<TH1F>> fiberLY_;

    std::vector<double> pixelsDataSummLY_[LAYERS_COUNT][CUBES_SLICE_X][CUBES_SLICE_Y];
    std::vector<double> pixelsDataXLY_[LAYERS_COUNT][CUBES_SLICE_X][CUBES_SLICE_Y];
    std::vector<double> pixelsDataYLY_[LAYERS_COUNT][CUBES_SLICE_X][CUBES_SLICE_Y];

    std::vector<std::unique_ptr<TH2F>> histogramsLY_gaus_layer_;
    std::vector<std::unique_ptr<TH2F>> x_histogramsLY_gaus_layer_;
    std::vector<std::unique_ptr<TH2F>> y_histogramsLY_gaus_layer_;

    std::vector<double> pixelsAverageSummLY_[CUBES_SLICE_X][CUBES_SLICE_Y];
    std::vector<double> pixelsAverageXLY_[CUBES_SLICE_X][CUBES_SLICE_Y];
    std::vector<double> pixelsAverageYLY_[CUBES_SLICE_X][CUBES_SLICE_Y];

    std::unique_ptr<TH2F> histogramsAverageLY_gaus_layer_;
    std::unique_ptr<TH2F> x_histogramsAverageLY_gaus_layer_;
    std::unique_ptr<TH2F> y_histogramsAverageLY_gaus_layer_;

    void NormalizeCubes() const;
    void NormalizeFibers() const;
    static void FillPixelData(TDirectory* dir, std::vector<double>* data, TF1* fit, TH2F* hist, const char* title, int x, int y, double cut);
    static double FillPixelData(std::vector<double>* data, TF1* fit, TH2F* hist, TH2F* histError, TH2F* histSigma, const char* title, int x, int y, double cut);
    static void FillPixelFiberData(TDirectory* dir, std::vector<double>* data, TH1F* hist, TH1F* histAvr, TH1F* histCompr, const char* title, int position);
    static double FillPixelProjectionXData(TDirectory* dirName, const std::vector<double> (&pixels)[LAYERS_COUNT][CUBES_SLICE_X][CUBES_SLICE_Y], TF1* fit, TH1F* hist, const char* title, int x_, int binX_min, int binX_max, int y_, int binY_min, int binY_max, int layer, double cut);
    static double FillPixelProjectionYData(TDirectory* dirName, const std::vector<double> (&pixels)[LAYERS_COUNT][CUBES_SLICE_X][CUBES_SLICE_Y], TF1* fit, TH1F* hist, const char* title, int x_, int binX_min, int binX_max, int y_, int binY_min, int binY_max, int layer, double cut);
};

#endif //CUBESLY_HISTOGRAMMANAGER_H
