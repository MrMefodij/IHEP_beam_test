//
// Created by amefodev on 05.11.2024.
//

#include <string>
#include <fstream>
#include "defines.h"
#include "connection_map.h"
#include "FileHandler.h"
#include "PrototypeHit.h"
#include "HistogramManager.h"

int main( int argc, char **argv ) {
    std::string filename = argv[1];
    Connection_Map* connectionMap;
        connectionMap = new Connection_Map();
        connectionMap->Init();
    const std::map<int,CubeConnection>* map = connectionMap->GetCubesMap();
    const std::map<int,CubePosition>* locationMap = connectionMap->GetCubesLocation();
    const std::map<int,GeomPosition>* channels_map = connectionMap->GetChannelMap();
    // const std::map<int, NormalizationConst>* cubes_normalization = connectionMap->GetCubeNormalization();
    const std::map<int, double>* fibers_normalization = connectionMap->GetFibersNormalization();

    LyHistogramManager hists;
    hists.SetTitlesHist("X[mm]","Y[mm]","N");
    hists.SetTitlesHistLY("X[mm]","Y[mm]","LY[p.e.]");
    hists.SetTitlesHistLY_Fiber_A("Amplitude LY[p.e.]","N");
    hists.SetTitlesHistLY_Fiber_I("Integral LY[p.e.]","N");
    hists.SetTitlesHistLY_Fiber_A_reconstructed("Reconstructed Amplitude LY[p.e.]","N");

    size_t pos = filename.rfind('/');
    std::string OutputDirName;
    if (pos != std::string::npos) {
        OutputDirName = filename.substr(0, pos+1);
    }
    std::shared_ptr<TFile> out = std::make_shared<TFile>((OutputDirName + "All_runs_LY.root").c_str(), "recreate");

    std::vector<std::string> dataFiles = GetDataFiles(OutputDirName, ".root", "_HitStructure.");
    for (const auto& file : dataFiles) {

        HitFileInput files(file, "cubes_position", false);
        PrototypeHit *hit = files.getHitPointer();

        for (int i = 0; i < files.GetEntries(); ++i) {
            files.GetNextEntry(i);
            for (int cube = 0; cube < CUBES_COUNT; ++cube) {
                auto itConnection = map->find(cube);
                if (itConnection == map->end()) continue;
                unsigned int x = itConnection->second._x;
                unsigned int y = itConnection->second._y;
                unsigned int z = itConnection->second._z;

                double z_position = channels_map->at(x)._cube_z;

                auto itLocation = locationMap->find(cube);
                if (itLocation == locationMap->end()) continue;
                double x_min = itLocation->second._x_min;
                double x_max = itLocation->second._x_max;
                double y_min = itLocation->second._y_min;
                double y_max = itLocation->second._y_max;


                auto xt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitX) {
                    return hitX._channel_id == x;
                });
                auto yt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitY) {
                    return hitY._channel_id == y;
                });
                auto zt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitZ) {
                    return hitZ._channel_id == z;
                });
                double x_position = hit->_Ax * (CUBES_POSITION - 1 + z_position)  + hit->_Bx;
                double y_position = hit->_Ay * (CUBES_POSITION - 1 + z_position) + hit->_By;

                if (xt != hit->_channelsData.end() && yt != hit->_channelsData.end() && zt != hit->_channelsData.end()) {
                    if (x_position >= x_min &&  x_position <= x_max+0.5 && y_position >= y_min && y_position <= y_max+0.5) {
                        if (xt->_chargePE_Amplitude > 1 && yt->_chargePE_Amplitude > 1 && zt->_amplitude > 1 && xt->_quality < 100 && yt->_quality < 100 ) {
                            // hists.fillHistograms(j,  x_position, y_position);
                            hists.SetBinContentLY(cube,  x_position, y_position, xt->_chargePE_Amplitude, yt->_chargePE_Amplitude);
                            // double x_ly = xt->_reconstructed_chargePE_Amplitude / cubes_normalization->at(cube).x_;
                            // double y_ly = yt->_reconstructed_chargePE_Amplitude / cubes_normalization->at(cube).y_;
                            double x_ly = xt->_chargePE_Amplitude / (fibers_normalization->at(xt->_channel_id) * CHARGE_FIT_CONSTANT);
                            double y_ly = yt->_chargePE_Amplitude / (fibers_normalization->at(yt->_channel_id) * CHARGE_FIT_CONSTANT);
                            if (xt->_amplitude < AMPLITUDE_OVERFLOW) {
                                x_ly = xt->_chargePE_Amplitude_XTalk / fibers_normalization->at(xt->_channel_id);
                            }
                            if (yt->_amplitude < AMPLITUDE_OVERFLOW) {
                                y_ly = yt->_chargePE_Amplitude_XTalk / fibers_normalization->at(yt->_channel_id);
                            }
                            hists.SetBinContentNormalizedLY(cube,  x_position, y_position, x_ly, y_ly);
                            hists.FillEventTime(xt->_time);
                            hists.FillEventTime(yt->_time);
                            hists.FillEventTime(zt->_time);

                            hists.SetFiberLY_A(cube, xt->_chargePE_Amplitude_XTalk, yt->_chargePE_Amplitude_XTalk);
                            hists.SetFiberLY_I(cube, xt->_chargePE_Integral_XTalk, yt->_chargePE_Integral_XTalk);
                            hists.SetFiberLY_A_recon(cube, xt->_chargePE_Amplitude, yt->_chargePE_Amplitude);
                            hists.FillFiberLY(xt->_channel_id, yt->_channel_id, x_ly, y_ly);
                            // hists.SetFiberLY_A_recon(cube,x_ly, y_ly);
                            hists.FillAveragePlots(x_position - x_min, y_position - y_min, x_ly, y_ly);
                            hists.SetBinContentNormalizedLY_layer(static_cast<unsigned int>(z_position),  x_position, y_position, x_ly, y_ly);
                        }
                    }
                    if (x_position >= locationMap->at(GLOBAL_CUBE)._x_min &&  x_position < locationMap->at(GLOBAL_CUBE)._x_max && y_position >= locationMap->at(GLOBAL_CUBE)._y_min && y_position < locationMap->at(GLOBAL_CUBE)._y_max) {
                        if (xt->_chargePE_Amplitude > 1 && yt->_chargePE_Amplitude > 1 && zt->_amplitude > 1 && xt->_quality < 100 && yt->_quality < 100) {
                            // double x_ly = xt->_reconstructed_chargePE_Amplitude / cubes_normalization->at(cube).x_;
                            // double y_ly = yt->_reconstructed_chargePE_Amplitude / cubes_normalization->at(cube).y_;
                            double x_ly = xt->_chargePE_Amplitude / (fibers_normalization->at(xt->_channel_id) * CHARGE_FIT_CONSTANT);
                            double y_ly = yt->_chargePE_Amplitude / (fibers_normalization->at(yt->_channel_id) * CHARGE_FIT_CONSTANT);
                            if (xt->_amplitude < AMPLITUDE_OVERFLOW) {
                                x_ly = xt->_chargePE_Amplitude_XTalk / fibers_normalization->at(xt->_channel_id);
                            }
                            if (yt->_amplitude < AMPLITUDE_OVERFLOW) {
                                y_ly = yt->_chargePE_Amplitude_XTalk / fibers_normalization->at(yt->_channel_id);
                            }
                            // hists.SetBinContentNormalizedLY_layer(static_cast<unsigned int>(z_position),  x_position, y_position, x_ly, y_ly);
                            hists.SetFiberLY_layer( xt->_channel_id, yt->_channel_id, x_position, y_position, x_ly, y_ly);
                        }
                    }
                }
            }
        }
    }
    hists.PreparePlotsForMC(out.get(), locationMap);
    hists.PreparePlotsForGeometry(out.get(), locationMap);
    hists.WriteHistograms(out.get());
    return 0;
}
