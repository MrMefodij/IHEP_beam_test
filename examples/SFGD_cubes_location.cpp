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

        // для работы с жёстким диском
        size_t lastSlash = file.rfind('/');
        std::string fileNameOnly = (lastSlash == std::string::npos) ? file : file.substr(lastSlash + 1);

        if (!fileNameOnly.empty() && fileNameOnly[0] == '.') {
            std::cout << "Пропускаем временный файл: " << file << std::endl;
            continue;
        }

        HitFileInput files(file, "cubes_position", false);
        PrototypeHit *hit = files.getHitPointer();

        for (int i = 0; i < files.GetEntries(); ++i) {
            files.GetNextEntry(i);
            if (!hit) continue;

            double layerLy[LAYERS_COUNT] = {0};
            for (const auto &itHit: hit->_channelsData) {
                if (channels_map->count(itHit._channel_id)) {
                    double norm_val = 1.0;
                    auto it_norm = fibers_s_normalization->find(itHit._channel_id);
                    if (it_norm != fibers_s_normalization->end()) norm_val = it_norm->second;
                    if (norm_val <= 0) norm_val = 1.0;

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
                double x_pos = hit->_Ax * (CUBES_POSITION - 1 + layer) + hit->_Bx;
                double y_pos = hit->_Ay * (CUBES_POSITION - 1 + layer) + hit->_By;

                // Попытка убрать точку (0,0), видимо, что-то не прошло инцилизацию
                if (hit->_Ax == 0 && hit->_Ay == 0) continue;

                hists.SetBinContentNormalizedLY_layer(layer, x_pos, y_pos, layerLy[layer], 0);
                hists.SetFiberLY_layer(layer, layer, x_pos, y_pos, layerLy[layer], layerLy[layer]);
            }

            for (int j = 0; j < CUBES_COUNT; ++j) {
                auto it_cube = map->find(j);
                if (it_cube == map->end()) continue;

                int x_id = it_cube->second._x;
                int y_id = it_cube->second._y;
                int z_id = it_cube->second._z;

                auto xt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(),
                                       [&](const auto &h) { return h._channel_id == x_id; });
                auto yt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(),
                                       [&](const auto &h) { return h._channel_id == y_id; });
                auto zt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(),
                                       [&](const auto &h) { return h._channel_id == z_id; });

                if (xt != hit->_channelsData.end() && yt != hit->_channelsData.end() &&
                    zt != hit->_channelsData.end()) {
//                    if (xt->_amplitude > 80 && yt->_amplitude > 80 && zt->_amplitude > 80) { // кат на порог в каналах ADC
                    if (xt->_chargePE_Amplitude > 8 && yt->_chargePE_Amplitude > 8 && zt->_chargePE_Amplitude > 8) { // кат на порог в PE
                        double z_pos_geom = channels_map->at(x_id)._cube_z;
                        double x_p = hit->_Ax * (CUBES_POSITION - 1 + z_pos_geom) + hit->_Bx;
                        double y_p = hit->_Ay * (CUBES_POSITION - 1 + z_pos_geom) + hit->_By;

                        if (hit->_Ax == 0 && hit->_Ay == 0) continue;

                        hists.fillHistograms(j, x_p, y_p);

                    }
                }
            }
        }
    }
//    hists.WriteHistograms(out.get());
    std::ofstream fileOut(CUBES_BOUNDARIES_MAP);
//    int threshold = 70;
//    hists.FindHighDensityArea(threshold, fileOut);
    hists.FindHighDensityArea(900, fileOut); // спец. кат на количество вхождений в гистограмму
    hists.WriteHistograms(out.get());
    return 0;
}
