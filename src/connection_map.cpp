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

Connection_Map::Connection_Map(std::string map_file, std::string cubes_map, std::string cubes_location, std::string channels_vs_cubes, std::string cubes_normalization, std::string fibers_normalization, std::string fibers_location, const unsigned int verbose)
: channelsFileName_(std::move(map_file)), cubesFileName_(std::move(cubes_map)), cubesLocationFileName_(std::move(cubes_location)),
channelsVsCubesFileName_(std::move(channels_vs_cubes)), cubes_normalizationFileName_(std::move(cubes_normalization)),
fibers_normalizationFileName_(std::move(fibers_normalization)), fibers_location_(std::move(fibers_location)),verbose_(verbose){}

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

    if(fibers_location_.find(".txt") != std::string::npos) FiberLocationByTXT();
    else {
        std::cerr << "Unexpected format of connection map" << std::endl;
    }
    std::cout << "Read: " << fibers_location_ << std::endl;
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

        unsigned int ch;
        try {
            ch = std::stoi(ch_str);
        } catch (...) {
            continue; // skip if channel is not a valid number
        }

        float x = std::nanf(""), y = std::nanf(""), z = std::nanf("");

        if (x_str != "NaN") x = std::stof(x_str);
        if (y_str != "NaN") y = std::stof(y_str);
        if (z_str != "NaN") z = std::stof(z_str);

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
    while (!map.eof()){
        unsigned int cube;
        map >> cube;
        if(std::floor(cube) == cube){
            map >> cubesMap_[cube]._x >> cubesMap_[cube]._y >> cubesMap_[cube]._z;
            if (verbose_){
                std::cout << AS_KV(cube) << " " << cubesMap_[cube]._x << " " << cubesMap_[cube]._y << " " <<  cubesMap_[cube]._z << std::endl;
            }
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
        unsigned int cube;
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

void Connection_Map::InitChannelsVsCubesByTXT(){
    std::ifstream map(channelsVsCubesFileName_.c_str());
    if(!map.is_open()){
        throw std::exception();
    }
    std::string header;
    std::getline(map, header);
    while (!map.eof()){
        unsigned int ch;
        map >> ch;
        if(std::floor(ch) == ch){
            map >> channelsVsCubes_[ch]._side >>  channelsVsCubes_[ch]._cube1 >>  channelsVsCubes_[ch]._cube2 >>  channelsVsCubes_[ch]._cube3;
            if (verbose_) {
                std::cout << AS_KV(ch) << " " << channelsVsCubes_[ch]._side << " " <<  channelsVsCubes_[ch]._cube1 << " " <<
                    channelsVsCubes_[ch]._cube2 << " " <<  channelsVsCubes_[ch]._cube3;
            }
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
        unsigned int cube;
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
        unsigned int cube;
        flag >> cube;
        if(std::floor(cube) == cube){
            flag >> fibersNormalization_[cube];
        }
    }
    fibersNormalization_Exist_ = true;
}

void Connection_Map::FiberLocationByTXT() {
    std::ifstream flag(fibers_location_.c_str());
    if(!flag.is_open()){
        std::cerr << fibers_location_ << " not found. Fibers location is not available.\n";
        return;
    }
    std::string header;
    std::getline(flag, header);
    while (!flag.eof()) {
        unsigned int ch;
        flag >> ch;
        if(std::floor(ch) == ch){
            flag >> fibersLocation_[ch].direction_ >> fibersLocation_[ch].x_min_ >> fibersLocation_[ch].x_max_ >> fibersLocation_[ch].y_min_ >> fibersLocation_[ch].y_max_;
        }
    }
}

const std::map<unsigned int,GeomPosition>* Connection_Map::GetChannelMap() const{
    return &channelsMap_;
}
const std::map<unsigned int,CubeConnection>* Connection_Map::GetCubesMap() const{
    return &cubesMap_;
}

const std::map<unsigned int,CubePosition>* Connection_Map::GetCubesLocation() const{
    return &cubesLocation_;
}

const std::map<unsigned int,ChannelsVsCubes>* Connection_Map::GetChannelsVsCubes() const{
    return &channelsVsCubes_;
}

const std::map<unsigned int, NormalizationConst>* Connection_Map::GetCubeNormalization() const{
    if (cubesNormalization_Exist_){
        return &cubesNormalization_;
    }
    return nullptr;
}

const std::map<unsigned int, double>* Connection_Map::GetFibersNormalization() const {
    if (fibersNormalization_Exist_){
        return &fibersNormalization_;
    }
    return nullptr;
}

const std::map<unsigned int, FibrePosition>* Connection_Map::GetFibersLocation() const {
    return &fibersLocation_;
}