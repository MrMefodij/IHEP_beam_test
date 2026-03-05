//
// Created by amefodev on 30.10.2024.
//

#ifndef CUBESLY_CONNECTION_MAP_H
#define CUBESLY_CONNECTION_MAP_H

#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include "defines.h"

struct GeomPosition{
    double _cube_x{NAN};
    double _cube_y{NAN};
    double _cube_z{NAN};
};

struct CubeConnection{
    int _x;
    int _y;
    int _z;
};

struct CubePosition{
    double _x_min{NAN};
    double _x_max{NAN};
    double _y_min{NAN};
    double _y_max{NAN};
};

struct ChannelsVsCubes {
    char _side;
    std::vector<int> cubes;
};

struct NormalizationConst {
    double x_;
    double y_;
    double z_;
};

struct FibrePosition {
    char direction_;
    double x_min_{NAN};
    double x_max_{NAN};
    double y_min_{NAN};
    double y_max_{NAN};
};

struct CubeMapLY {
    int x_;
    int y_;
    double ly_;
    bool operator==(const CubeMapLY& other) const {
        return x_ == other.x_ && y_ == other.y_;
    }
};

class Connection_Map {
public:
    explicit Connection_Map(std::string map_file = CHANNELS_MAP, std::string cubes_map = CUBES_MAP, std::string cubes_location = CUBES_BOUNDARIES_MAP,
                            std::string channels_vs_cubes = CHANNELS_VS_CUBES, std::string cubes_normalization = CUBES_NORMALIZATION,
                            std::string fibers_normalization = FIBERS_NORMALIZATION, int verbose=0);
    void Init();
    const std::map<int,GeomPosition>* GetChannelMap() const;
    const std::map<int,CubeConnection>* GetCubesMap() const;
    const std::map<int,CubePosition>* GetCubesLocation() const;
    const std::map<int, ChannelsVsCubes>* GetChannelsVsCubes() const;
    const std::map<int, NormalizationConst>* GetCubeNormalization() const;
    const std::map<int, double>* GetFibersNormalization() const;
    const std::map<int, FibrePosition>* GetFibersLocation() const;

private:
    void InitGeomByTXT();
    void InitCubesByTXT();
    void InitLocationByTXT();
    void InitChannelsVsCubesByTXT();
    void InitCubesNormalizationByTXT();
    void InitFibersNormalizationByTXT();
    void FiberLocationByTXT();

    const std::string channelsFileName_;
    const std::string cubesFileName_;
    const std::string cubesLocationFileName_;
    const std::string channelsVsCubesFileName_;
    const std::string cubes_normalizationFileName_;
    const std::string fibers_normalizationFileName_;
    int verbose_;
    std::map<int,GeomPosition> channelsMap_;
    std::map<int,CubeConnection> cubesMap_;
    std::map<int,CubePosition> cubesLocation_;
    std::map<int,ChannelsVsCubes> channelsVsCubes_;

    bool cubesNormalization_Exist_;
    std::map<int, NormalizationConst> cubesNormalization_;

    bool fibersNormalization_Exist_;
    std::map<int, double> fibersNormalization_;
};


#endif //CUBESLY_CONNECTION_MAP_H
