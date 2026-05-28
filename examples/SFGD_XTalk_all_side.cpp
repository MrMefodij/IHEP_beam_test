//
// Created by Angelina Chvirova on 23.01.2025.
//
#include <string>
#include <fstream>
#include <memory>
#include <TFile.h>
#include "defines.h"
#include "connection_map.h"
#include "FileHandler.h"
#include "PrototypeHit.h"
#include "HistogramManager.h"

template<typename T>
bool Find_cube(const std::map<int,CubeConnection>* map, const std::map<int,CubePosition>* locationMap, const std::map<int,GeomPosition>* channels_map,
               int& cube,  PrototypeHit* hit, T& xt, T& yt, T& zt, double& x_min, double& x_max, double& y_min, double& y_max, double& x_position_cube, double& y_position_cube,  double& z_position_cube){
    auto itConnection = map->find(cube);
    if (itConnection == map->end()) return false;
    int x = itConnection->second._x;
    int y = itConnection->second._y;
    int z = itConnection->second._z;

    // если канал не существует (-1)
    if (x < 0 || y < 0 || z < 0) return false;

    if (channels_map->find(x) == channels_map->end()) {
        std::cerr << "Warning: Channel " << x << " not in geometry map!" << std::endl;
        return false;
    }

    x_position_cube = channels_map->at(y)._cube_x;
    y_position_cube = channels_map->at(x)._cube_y;
    z_position_cube = channels_map->at(x)._cube_z;

    auto itLocation = locationMap->find(cube);
    if (itLocation == locationMap->end()) return false;
    x_min = itLocation->second._x_min;
    x_max = itLocation->second._x_max;
    y_min = itLocation->second._y_min;
    y_max = itLocation->second._y_max;


    xt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitX) {
        return hitX._channel_id == x;
    });
    yt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitY) {
        return hitY._channel_id == y;
    });
    zt = std::find_if(hit->_channelsData.begin(), hit->_channelsData.end(), [&](const auto &hitZ) {
        return hitZ._channel_id == z;
    });
    return true;
}

int main( int argc, char **argv ) {
    std::string filename = argv[1];

    Connection_Map connectionMap;
    connectionMap.Init();
    const std::map<int,CubeConnection>* map = connectionMap.GetCubesMap();
    const std::map<int,CubePosition>* locationMap = connectionMap.GetCubesLocation();
    const std::map<int,GeomPosition>* channels_map = connectionMap.GetChannelMap();
    const std::map<int, double>* fibers_s_normalization = connectionMap.GetFibersNormalization();

    XTalk4sideHistogramManager hists;
    hists.SetTitlesHistXTalk4side_1D("CrossTalk", "N");
    hists.SetTitlesHistXTalk4side_1D_2("X or Y [mm]", "CrossTalk");
    hists.SetTitlesHistXTalk4side_1D_3("LY [p.e.]", "N");
    hists.SetTitlesHistXTalk4side_2D("X [mm]","Y [mm]", "CrossTalk");

    size_t pos = filename.rfind('/');
    std::string OutputDirName;
    if (pos != std::string::npos) {
        OutputDirName = filename.substr(0, pos + 1);
    }
    std::shared_ptr<TFile> out = std::make_shared<TFile>((OutputDirName + "All_XTalk_4_side.root").c_str(), "recreate");

    std::vector<std::string> dataFiles = GetDataFiles(OutputDirName, ".root", "_HitStructure.");
    for (const auto& file : dataFiles) {

        HitFileInput files(file, "cubes_position", false);
        PrototypeHit *hit = files.getHitPointer();

        for (int i = 0; i < files.GetEntries(); ++i) {
            files.GetNextEntry(i);
            for (int cube = 0; cube < CUBES_COUNT; ++cube) {
                auto xt = hit->_channelsData.begin();
                auto yt = hit->_channelsData.begin();
                auto zt = hit->_channelsData.begin();
                double x_min, x_max, y_min, y_max, x_position_cube, y_position_cube, z_position_cube;
                bool found = Find_cube(map, locationMap, channels_map, cube, hit, xt, yt, zt, x_min, x_max, y_min,
                                       y_max, x_position_cube, y_position_cube, z_position_cube);
                if (!found) continue;

                // сначала нахожу куб, где прошла частица
                double x_position = hit->_Ax * (CUBES_POSITION - 1 + z_position_cube) + hit->_Bx;
                double y_position = hit->_Ay * (CUBES_POSITION - 1 + z_position_cube) + hit->_By;

                // частица ГЕОМЕТРИЧЕСКИ попал в текущий куб
                if (x_position >= x_min && x_position < x_max && y_position >= y_min && y_position < y_max) {

                    if (xt != hit->_channelsData.end() && yt != hit->_channelsData.end() &&
                        zt != hit->_channelsData.end()) {
                        auto it_x = fibers_s_normalization->find(xt->_channel_id);
                        auto it_y = fibers_s_normalization->find(yt->_channel_id);

                        if (it_x == fibers_s_normalization->end() || it_y == fibers_s_normalization->end()) {
                            std::cerr << "Ошибка: канал не найден в fibers_s_normalization" << std::endl;
                            std::cerr << "Файл: " << file << std::endl;
                            if (it_x == fibers_s_normalization->end()) {
                                std::cerr << "Не найден канал с id: " << xt->_channel_id << std::endl;
                            }
                            if (it_y == fibers_s_normalization->end()) {
                                std::cerr << "Не найден канал с id: " << yt->_channel_id << std::endl;
                            }
                            continue;
                        }

                        double x_ly = xt->_chargePE_Amplitude / fibers_s_normalization->at(xt->_channel_id);
                        double y_ly = yt->_chargePE_Amplitude / fibers_s_normalization->at(yt->_channel_id);

                        if (xt != hit->_channelsData.end() && yt != hit->_channelsData.end() &&
                            zt != hit->_channelsData.end()) {
                            if ((x_ly > 0 || y_ly > 0)
                                && xt->_quality < 100 && yt->_quality < 100 && xt->_time > 205 && xt->_time < 280 &&
                                yt->_time > 205 && yt->_time < 280 && zt->_time > 205 && zt->_time < 280) {

                                for (int cube_xTalk: {cube - 21, cube + 21, cube - 3, cube + 3}) {

                                    if (cube_xTalk < 1 || cube_xTalk > 102) continue; // 102 - зависит от размеров сборки

                                    auto itConnection = map->find(cube_xTalk);
                                    if (itConnection == map->end()) continue;

                                    // каналы
                                    int x_xTalk = itConnection->second._x;
                                    int y_xTalk = itConnection->second._y;
                                    int z_xTalk = itConnection->second._z;

                                    // если канал не существует (-1)
                                    if (x_xTalk < 0 || y_xTalk < 0 || z_xTalk < 0) continue;

                                    auto xt_xTalk = std::find_if(hit->_channelsData.begin(),
                                                                 hit->_channelsData.end(),
                                                                 [&](const auto &ch) {
                                                                     return ch._channel_id == x_xTalk;
                                                                 });
                                    auto yt_xTalk = std::find_if(hit->_channelsData.begin(),
                                                                 hit->_channelsData.end(),
                                                                 [&](const auto &ch) {
                                                                     return ch._channel_id == y_xTalk;
                                                                 });

                                    auto zt_xTalk = std::find_if(hit->_channelsData.begin(),
                                                                 hit->_channelsData.end(),
                                                                 [&](const auto &ch) {
                                                                     return ch._channel_id == z_xTalk;
                                                                 });

                                    double x_ly_xTalk = 0.0;
                                    double y_ly_xTalk = 0.0;

                                    if (xt_xTalk->_channel_id == -1 || yt_xTalk->_channel_id == -1 || zt_xTalk->_channel_id == -1) continue;

                                    if (xt_xTalk != hit->_channelsData.end() && (xt_xTalk->_time > 205 && xt_xTalk->_time < 280)){
                                        x_ly_xTalk = xt_xTalk->_chargePE_Amplitude / fibers_s_normalization->at(xt_xTalk->_channel_id);
                                    }

                                    if (yt_xTalk != hit->_channelsData.end() && (yt_xTalk->_time > 205 && yt_xTalk->_time < 280)) {
                                        y_ly_xTalk = yt_xTalk->_chargePE_Amplitude / fibers_s_normalization->at(yt_xTalk->_channel_id);
                                    }

//                                    if (x_ly_xTalk >= 0 && y_ly_xTalk >= 0) {
                                        if ((cube_xTalk - cube) == -21) {
                                            hists.SetBinContentXTalk4side_up(x_position - x_min,
                                                                             y_position - y_min,
                                                                             x_ly_xTalk / x_ly, x_ly,
                                                                             x_ly_xTalk);
                                        }
                                        if ((cube_xTalk - cube) == 21) {
                                            hists.SetBinContentXTalk4side_down(x_position - x_min,
                                                                               y_position - y_min,
                                                                               x_ly_xTalk / x_ly, x_ly,
                                                                               x_ly_xTalk);
                                        }

                                        if ((cube_xTalk - cube) == 3) {
                                            hists.SetBinContentXTalk4side_left(x_position - x_min,
                                                                               y_position - y_min,
                                                                               y_ly_xTalk / y_ly, y_ly,
                                                                               y_ly_xTalk);
                                        }
                                        if ((cube_xTalk - cube) == -3) {
                                            hists.SetBinContentXTalk4side_right(x_position - x_min,
                                                                                y_position - y_min,
                                                                                y_ly_xTalk / y_ly, y_ly,
                                                                                y_ly_xTalk);
                                        }
//                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    hists.WriteHistograms(out.get());
    return 0;
}