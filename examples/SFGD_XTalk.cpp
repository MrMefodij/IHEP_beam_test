//
// Created by amefodev on 07.12.2024.
//

#include <string>
#include <fstream>
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

    // Если канал не существует (-1)
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
    const std::map<int, double> *fibers_normalization = connectionMap.GetFibersNormalization();

    XTalkHistogramManager hists;
    hists.SetTitlesHistXTalk_N("X or Y[mm]","N");
    hists.SetTitlesHistXTalk_1D("X or Y[mm]","CrossTalk");
    hists.SetTitlesHistXTalk_1D_all("CrossTalk", "N");
    hists.SetTitlesHistXTalk_2D("X[mm]","Y[mm]", "CrossTalk");

    size_t pos = filename.rfind('/');
    std::string OutputDirName;
    if (pos != std::string::npos) {
        OutputDirName = filename.substr(0, pos + 1);
    }
    std::shared_ptr<TFile> out = std::make_shared<TFile>((OutputDirName + "All_XTalk.root").c_str(), "recreate");

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

                double x_position = hit->_Ax * (CUBES_POSITION - 1 + z_position_cube) + hit->_Bx;
                double y_position = hit->_Ay * (CUBES_POSITION - 1 + z_position_cube) + hit->_By;

//                if (x_ly < 0 || y_ly < 0) {
//                    std::cout << x_ly << " " << y_ly << std::endl;
//                }

                if (xt != hit->_channelsData.end() && yt != hit->_channelsData.end() &&
                    zt != hit->_channelsData.end()) {

                    double x_ly = xt->_chargePE_Amplitude / (fibers_normalization->at(xt->_channel_id));
                    double y_ly = yt->_chargePE_Amplitude / (fibers_normalization->at(yt->_channel_id));

                    if (fibers_normalization->find(cube) == fibers_normalization->end()) continue;

                    if (xt->_chargePE_Amplitude >= 0 && yt->_chargePE_Amplitude >= 0 &&
                        zt->_amplitude > 1 // обсудить: нужен ли этот кат?
                        && xt->_quality < 100 && yt->_quality < 100 && xt->_time > 205 && xt->_time < 280 &&
                        yt->_time > 205 && yt->_time < 280 && zt->_time > 205 && zt->_time < 280) {

                        for (int cube_xTalk = 0; cube_xTalk < CUBES_COUNT; ++cube_xTalk) {
                            auto xt_xTalk = hit->_channelsData.begin();
                            auto yt_xTalk = hit->_channelsData.begin();
                            auto zt_xTalk = hit->_channelsData.begin();
                            double x_min_xTalk, x_max_xTalk, y_min_xTalk, y_max_xTalk, x_position_xTalk, y_position_xTalk, z_position_xTalk;
                            bool found_xTalk = Find_cube(map, locationMap, channels_map, cube_xTalk, hit, xt_xTalk,
                                                         yt_xTalk, zt_xTalk, x_min_xTalk, x_max_xTalk, y_min_xTalk,
                                                         y_max_xTalk, x_position_xTalk, y_position_xTalk,
                                                         z_position_xTalk);
                            if (!found_xTalk) continue;

                            if (z_position_cube == z_position_xTalk && cube_xTalk != cube &&
                                xt_xTalk->_chargePE_Amplitude > 0 && yt_xTalk->_chargePE_Amplitude > 0 &&
                                xt_xTalk != hit->_channelsData.end() && yt_xTalk != hit->_channelsData.end() &&
                                zt_xTalk != hit->_channelsData.end()
                                && x_position >= x_min_xTalk && x_position < x_max_xTalk && y_position >= y_min_xTalk &&
                                y_position < y_max_xTalk
                                && xt_xTalk->_time > 205 && xt_xTalk->_time < 280 && yt_xTalk->_time > 205 &&
                                yt_xTalk->_time < 280 && zt_xTalk->_time > 205 && zt_xTalk->_time < 280) {

                                double x_ly_xTalk = xt_xTalk->_chargePE_Amplitude / (fibers_normalization->at(xt->_channel_id));
                                double y_ly_xTalk = yt_xTalk->_chargePE_Amplitude / (fibers_normalization->at(yt->_channel_id));

                                if (!std::isnan(x_ly) && !std::isnan(x_ly_xTalk) && !std::isnan(y_ly) &&!std::isnan(y_ly_xTalk)) {
                                    if (!(y_position >= y_min && y_position < y_max) &&
                                        x_position_cube == x_position_xTalk && (x_ly == 0 || x_ly > 0) && x_ly_xTalk > 15) {
                                        //                                  std::cout << cube << " " << " " << x_position_cube  << " " << y_position_cube <<  " " <<  z_position_cube << " " << cube_xTalk <<  " " << x_position_xTalk  << " " << y_position_xTalk <<  " " <<  z_position_xTalk << std::endl;
                                        hists.SetBinContentXTalk_y(cube, x_position, y_position,
                                                                   y_ly / y_ly_xTalk, x_ly / x_ly_xTalk);
                                        }

                                    if (!(x_position >= x_min && x_position < x_max) &&
                                        y_position_cube == y_position_xTalk && (y_ly == 0 || y_ly > 0) &&  y_ly_xTalk > 15) {
                                        //                                  std::cout << cube << " " << " " << x_position_cube  << " " << y_position_cube <<  " " <<  z_position_cube << " " << cube_xTalk <<  " " << x_position_xTalk  << " " << y_position_xTalk <<  " " <<  z_position_xTalk << std::endl;
                                        hists.SetBinContentXTalk_x(cube, x_position, y_position, y_ly / y_ly_xTalk,
                                                                   x_ly / x_ly_xTalk);
                                        }
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
