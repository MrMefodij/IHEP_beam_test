//
// Created by amefodev on 14.02.2025.
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
    bool angle = false;
    if (filename.find("testbench_61") != std::string::npos) {
        std::cout << "Applying angle cubes location\n";
        angle = true;
        connectionMap = new Connection_Map(CHANNELS_MAP_ANGLE, CUBES_MAP,CUBES_BOUNDARIES_MAP_ANGLE,CHANNELS_VS_CUBES,CUBES_NORMALIZATION_ANGLE, FIBERS_NORMALIZATION, FIBERS_BOUNDARIES_MAP, 0);
    } else {
        connectionMap = new Connection_Map();
    }
    connectionMap->Init();
    const std::map<unsigned int,CubeConnection>* map = connectionMap->GetCubesMap();
    const std::map<unsigned int,GeomPosition>* channels_map = connectionMap->GetChannelMap();
    const std::map<unsigned int, double>* fibers_s_normalization = connectionMap->GetFibersNormalization();
    LyHistogramManager hists;
    hists.SetTitlesHist("X[mm]","Y[mm]","Number");

    size_t pos = filename.rfind('/');
    std::string OutputDirName;
    if (pos != std::string::npos) {
        OutputDirName = filename.substr(0, pos+1);
    }

    std::shared_ptr<TFile> out = std::make_shared<TFile>((OutputDirName + "All_runs_cubes_position.root").c_str(), "recreate");

    std::vector<std::string> dataFiles = GetDataFiles(OutputDirName, ".root", "_HitStructure.");
    for (const auto& file : dataFiles) {

        HitFileInput files(file, "cubes_position", false);
        PrototypeHit *hit = files.getHitPointer();

        for (int i = 0; i < files.GetEntries(); ++i) {
            files.GetNextEntry(i);
            for (int j = 0; j < CUBES_COUNT; ++j) {
                double layerLy[3] = {0,0,0};
                double layerQuality[3] = {0,0,0};
                for (const auto & itHit : hit->_channelsData) {
                    if (channels_map->find(itHit._channel_id) != channels_map->end() &&  fibers_s_normalization->find(itHit._channel_id) != fibers_s_normalization->end()) {
                        auto ly = itHit._reconstructed_chargePE_Amplitude / fibers_s_normalization->at(itHit._channel_id);
                        if (!std::isnan(channels_map->at(itHit._channel_id)._cube_z)) {
                            layerLy[static_cast<int>(channels_map->at(itHit._channel_id)._cube_z)] += ly;
                        }
                    }
                }
                for (int layer = 0; layer < 3; ++layer) {
                    double x_position = hit->_Ax * (CUBES_POSITION - 1 + layer)  + hit->_Bx;
                    double y_position = hit->_Ay * (CUBES_POSITION - 1 + layer) + hit->_By;
                    if (layerQuality[layer] < 100) {
                        hists.SetBinContentNormalizedLY_layer(layer,  x_position, y_position, layerLy[layer], 0);
                        hists.SetFiberLY_layer(layer, layer, x_position, y_position, layerLy[layer], layerLy[layer]);
                    }
                }

                auto it = map->find(j);
                if (it == map->end()) continue;
                unsigned int x = it->second._x;
                unsigned int y = it->second._y;
                unsigned int z = it->second._z;

                double z_position = channels_map->at(x)._cube_z;
                if (angle) {
                    z_position = channels_map->at(j)._cube_z;
                }

                auto xt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitX) {
                    return hitX._channel_id == x;
                });
                auto yt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitY) {
                    return hitY._channel_id == y;
                });
                auto zt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitZ) {
                    return hitZ._channel_id == z;
                });
                if (xt != hit->_channelsData.end() && yt != hit->_channelsData.end() &&
                    zt != hit->_channelsData.end()) {
                    if (xt->_reconstructed_chargePE_Amplitude > 5 && yt->_reconstructed_chargePE_Amplitude > 5 && zt->_chargeADC > 1000 ) {
                        double x_position = hit->_Ax * (CUBES_POSITION - 1 + z_position)  + hit->_Bx;
                        double y_position = hit->_Ay * (CUBES_POSITION - 1 + z_position) + hit->_By;
                        if (angle) {
                            x_position = hit->_Ax * (CUBES_POSITION + z_position)  + hit->_Bx;
                            y_position = hit->_Ay * (CUBES_POSITION + z_position) + hit->_By;
                        }
                        hists.fillHistograms(j, x_position, y_position);
                    }
                }
            }
        }
    }
    hists.WriteHistograms(out.get());
    std::ofstream fileOut(CUBES_BOUNDARIES_MAP);
    int threshold = 80;
    if (angle) {
        threshold = 120;
    }
    hists.FindHighDensityArea(threshold, fileOut);
    return 0;
}
