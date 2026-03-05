//
// Created by amefodev on 30.10.2024.
//

#include "connection_map.h"
#include <exception>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include "PrototypeHit.h"

ClassImp(ChannelData)
ClassImp(PrototypeHit)

Connection_Map::Connection_Map(std::string map_file, std::string cubes_map, std::string cubes_location, std::string channels_vs_cubes, std::string cubes_normalization, std::string fibers_normalization, const int verbose)
: channelsFileName_(std::move(map_file)), cubesFileName_(std::move(cubes_map)), cubesLocationFileName_(std::move(cubes_location)),
channelsVsCubesFileName_(std::move(channels_vs_cubes)), cubes_normalizationFileName_(std::move(cubes_normalization)),
fibers_normalizationFileName_(std::move(fibers_normalization)), verbose_(verbose){}

void Connection_Map::Init(){
    if(channelsFileName_.find(".txt") != std::string::npos) InitGeomByTXT();
    else {
        std::cerr << "Unexpected format of connection map" << std::endl;
    }
    std::cout << "Read: " << channelsFileName_ << std::endl;

    if(cubesFileName_.find(".txt") != std::string::npos) InitCubesByTXT();
    else {
        std::cerr << "Unexpected format of connection map" << std::endl;
    }
    std::cout << "Read: " << cubesFileName_ << std::endl;

    if(cubesLocationFileName_.find(".txt") != std::string::npos) InitLocationByTXT();
    else {
        std::cerr << "Unexpected format of connection map" << std::endl;
    }
    std::cout << "Read: " << cubesLocationFileName_ << std::endl;

     if(channelsVsCubesFileName_.find(".txt") != std::string::npos) InitChannelsVsCubesByTXT();
     else {
         std::cerr << "Unexpected format of connection map" << std::endl;
     }
        std::cout << "Read: " << channelsVsCubesFileName_ << std::endl;

    if(cubes_normalizationFileName_.find(".txt") != std::string::npos) InitCubesNormalizationByTXT();
    else {
        std::cerr << "Unexpected format of connection map" << std::endl;
    }
    std::cout << "Read: " << cubes_normalizationFileName_ << std::endl;

    if(fibers_normalizationFileName_.find(".txt") != std::string::npos) InitFibersNormalizationByTXT();
    else {
        std::cerr << "Unexpected format of connection map" << std::endl;
    }
    std::cout << "Read: " << fibers_normalizationFileName_ << std::endl;
}

void Connection_Map::InitGeomByTXT(){
    std::ifstream map(channelsFileName_);
    if (!map.is_open()) {
        throw std::runtime_error("Cannot open file " + channelsFileName_);
    }

    std::string line;
    std::getline(map, line); // skip header

    while (std::getline(map, line)) {
        std::istringstream iss(line);
        std::string ch_str, x_str, y_str, z_str;
        if (!(iss >> ch_str >> x_str >> y_str >> z_str)) {
            continue; // skip malformed lines
        }

        int ch;
        try {
            ch = std::stoi(ch_str);
        } catch (...) {
            continue; // skip if channel is not a valid number
        }

        float x = -1.0, y = -1.0, z = -1.0;

        if (x_str != "-1") x = std::stof(x_str);
        if (y_str != "-1") y = std::stof(y_str);
        if (z_str != "-1") z = std::stof(z_str);

        channelsMap_[ch]._cube_x = x;
        channelsMap_[ch]._cube_y = y;
        channelsMap_[ch]._cube_z = z;

        if (verbose_) {
            std::cout << "ch: " << ch << " " << x << " " << y << " " << z << std::endl;
        }
    }
}

void Connection_Map::InitCubesByTXT(){
    std::ifstream map(cubesFileName_.c_str());
    if(!map.is_open()){
        throw std::exception();
    }
    std::string header;
    std::getline(map, header);
        int cube;

    while (map >> cube) {
        map >> cubesMap_[cube]._x;
        map >> cubesMap_[cube]._y;
        map >> cubesMap_[cube]._z;

        if (verbose_) {
            std::cout << "Read cube: " << cube << std::endl;
        }
    }
}

void Connection_Map::InitLocationByTXT(){
    std::ifstream map(cubesLocationFileName_.c_str());
    if(!map.is_open()){
        throw std::exception();
    }
    std::string header;
    std::getline(map, header);
    while (!map.eof()){
        int cube;
        map >> cube;
        if(std::floor(cube) == cube){
            map >> cubesLocation_[cube]._x_min >> cubesLocation_[cube]._x_max >> cubesLocation_[cube]._y_min >> cubesLocation_[cube]._y_max;
            if (verbose_){
                std::cout << AS_KV(cube) << " " << cubesLocation_[cube]._x_min << " " << cubesLocation_[cube]._x_max <<
                " " << cubesLocation_[cube]._y_min << " " << cubesLocation_[cube]._y_max << std::endl;
            }
        }
    }
}

void Connection_Map::InitChannelsVsCubesByTXT() {
    std::ifstream mapFile(channelsVsCubesFileName_.c_str());
    if (!mapFile.is_open()) throw std::runtime_error("File not found");

    std::string line;
    std::getline(mapFile, line);

    while (std::getline(mapFile, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        int ch;
        char side;

        if (!(ss >> ch >> side)) continue;

        channelsVsCubes_[ch]._side = side;
        channelsVsCubes_[ch].cubes.clear();

        int temp_cube;
        while (ss >> temp_cube) {
            channelsVsCubes_[ch].cubes.push_back(temp_cube);
        }

        if (verbose_) {
            std::cout << "Channel " << ch << " (" << side << ") has "
                      << channelsVsCubes_[ch].cubes.size() << " cubes." << std::endl;
        }
    }
}

void Connection_Map::InitCubesNormalizationByTXT(){
    std::ifstream flag(cubes_normalizationFileName_.c_str());
    if(!flag.is_open()){
        std::cerr << cubes_normalizationFileName_ << " not found. Cube flag is not available.\n";
        cubesNormalization_Exist_ = false;
        return;
    }
    for (int i = 0; i < CUBES_COUNT; ++i) {
        cubesNormalization_[i] = {1,1,1};
    }
    std::string header;
    std::getline(flag, header);
    while (!flag.eof()) {
        int cube;
        flag >> cube;
        if(std::floor(cube) == cube) {
            char side;
            flag >> side;
            if (side == 'x'){
                flag >> cubesNormalization_[cube].x_;
            } else if (side == 'y'){
                flag >> cubesNormalization_[cube].y_;
            } else if (side == 'z'){
                flag >> cubesNormalization_[cube].z_;
            } else{
                cubesNormalization_Exist_ = false;
                return;
            }
        }
    }
    cubesNormalization_Exist_ = true;
}

void Connection_Map::InitFibersNormalizationByTXT() {
    std::ifstream flag(fibers_normalizationFileName_.c_str());
    if(!flag.is_open()){
        std::cerr << fibers_normalizationFileName_ << " not found. Channels flag is not available.\n";
        fibersNormalization_Exist_ = false;
        return;
    }
    std::string header;
    std::getline(flag, header);
    while (!flag.eof()) {
        int cube;
        flag >> cube;
        if(std::floor(cube) == cube){
            flag >> fibersNormalization_[cube];
        }
    }
    fibersNormalization_Exist_ = true;
}

const std::map<int,GeomPosition>* Connection_Map::GetChannelMap() const{
    return &channelsMap_;
}
const std::map<int,CubeConnection>* Connection_Map::GetCubesMap() const{
    return &cubesMap_;
}

const std::map<int,CubePosition>* Connection_Map::GetCubesLocation() const{
    return &cubesLocation_;
}

const std::map<int,ChannelsVsCubes>* Connection_Map::GetChannelsVsCubes() const{
    return &channelsVsCubes_;
}

const std::map<int, NormalizationConst>* Connection_Map::GetCubeNormalization() const{
    if (cubesNormalization_Exist_){
        return &cubesNormalization_;
    }
    return nullptr;
}

const std::map<int, double>* Connection_Map::GetFibersNormalization() const {
    if (fibersNormalization_Exist_){
        return &fibersNormalization_;
    }
    return nullptr;
}