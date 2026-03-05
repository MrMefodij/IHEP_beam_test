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
    std::cout << "Starting Init..." << std::endl;
    connectionMap = new Connection_Map();
    connectionMap->Init();
    std::cout << "Init Finished!" << std::endl;
    const std::map<int,CubeConnection>* map = connectionMap->GetCubesMap();
    const std::map<int,GeomPosition>* channels_map = connectionMap->GetChannelMap();
    const std::map<int, double>* fibers_s_normalization = connectionMap->GetFibersNormalization();
    LyHistogramManager hists;
    hists.SetTitlesHist("X[mm]","Y[mm]","Number");

    size_t pos = filename.rfind('/');
    std::string OutputDirName;
    if (pos != std::string::npos) {
        OutputDirName = filename.substr(0, pos+1);
    }

    std::shared_ptr<TFile> out = std::make_shared<TFile>((OutputDirName + "All_runs_cubes_position.root").c_str(), "recreate");

    std::vector<std::string> dataFiles = GetDataFiles(OutputDirName, ".root", "_HitStructure.");
    std::cout << "Files found: " << dataFiles.size() << std::endl;

    for (const auto& file : dataFiles) {

        HitFileInput files(file, "cubes_position", false);
        PrototypeHit *hit = files.getHitPointer();

        for (int i = 0; i < files.GetEntries(); ++i) {
            if (!channels_map || !fibers_s_normalization) {
                std::cerr << "Error: Nullptr!" << std::endl;
                return 1;
            }
            files.GetNextEntry(i);

            if (!hit) continue;

            for (int j = 0; j < CUBES_COUNT; ++j) {
                double layerLy[LAYERS_COUNT] = {0};
                double layerQuality[LAYERS_COUNT] = {0};

                for (const auto & itHit : hit->_channelsData) {
                    if (channels_map->count(itHit._channel_id)) {

                        double norm_val = 1.0;
                        if (fibers_s_normalization && fibers_s_normalization->count(itHit._channel_id)) {
                            norm_val = fibers_s_normalization->at(itHit._channel_id);
                        }

                        if (norm_val == 0) norm_val = 1.0;

                        auto ly = itHit._chargePE_Amplitude / norm_val;
                        double z_coord = channels_map->at(itHit._channel_id)._cube_z;

                        if (!std::isnan(z_coord) && z_coord >= 0 && z_coord < LAYERS_COUNT) {
                            int layer_idx = static_cast<int>(z_coord);
                            layerLy[layer_idx] += ly;
                        }
                    }
                }

                for (int layer = 0; layer < LAYERS_COUNT; ++layer) {
                    double x_position = hit->_Ax * (CUBES_POSITION - 1 + layer)  + hit->_Bx;
                    double y_position = hit->_Ay * (CUBES_POSITION - 1 + layer) + hit->_By;
                    if (layerQuality[layer] < 100) {
                        hists.SetBinContentNormalizedLY_layer(layer,  x_position, y_position, layerLy[layer], 0);
                        hists.SetFiberLY_layer(layer, layer, x_position, y_position, layerLy[layer], layerLy[layer]);
                    }
                }

                auto it = map->find(j);
                if (it == map->end()) continue;
                int x = static_cast<int>(it->second._x);
                int y = static_cast<int>(it->second._y);
                int z = static_cast<int>(it->second._z);

                if (x < 0 || y < 0 || z < 0) continue;
                if (channels_map->find(x) == channels_map->end() ||
                    channels_map->find(y) == channels_map->end() ||
                    channels_map->find(z) == channels_map->end()) {
                    continue;
                }

                double z_position = channels_map->at(x)._cube_z;

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
                    if (xt->_amplitude > 5 && yt->_amplitude > 5 && zt->_amplitude > 15) {
                        double x_position = hit->_Ax * (CUBES_POSITION - 1 + z_position)  + hit->_Bx;
                        double y_position = hit->_Ay * (CUBES_POSITION - 1 + z_position) + hit->_By;
                        hists.fillHistograms(j, x_position, y_position);
                    }
                }
            }
        }
    }
    hists.WriteHistograms(out.get());
    std::ofstream fileOut(CUBES_BOUNDARIES_MAP);
    int threshold = 70;
    hists.FindHighDensityArea(threshold, fileOut);
    return 0;
}
