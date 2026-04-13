//
// Created by amefodev on 30.10.2024.
//

#ifndef CUBESLY_PNPI_STRUCTURE_H
#define CUBESLY_PNPI_STRUCTURE_H

#include <iostream>
#include <cmath>
#include <set>
#include "defines.h"
#include "TFile.h"
#include "TTree.h"
#include "PrototypeHit.h"
#include "Calibration.h"

namespace PNPI {
    struct WaveFormParamStruct {
        float A;    // amplitude
        float T;    // time
        float I;    // integral
        float Q;    // quality
        float B;    // baseline
        float AF;   //reconstructed A амплитуда из фита, если фит сошелся. если не сошелся -- отрицательное число
        int overflow; //было ли переполнение по амплитуде
    };

    struct TrackParamStruct {
        float AX;    // track X: a*x + b
        float BX;
        float AY;    // track Y: a*y + b
        float BY;
        float Q;    // quality - negative if bad track in one of the projections
    };

    struct TreeStructure{
        TrackParamStruct _trackParam;
        std::map<int,WaveFormParamStruct> _sipmData;
    };

    PrototypeHit GetHitStructure(const TreeStructure* entry, const std::set<int>* availableChannels_, const Calibration* calibration = nullptr);

    class PnpiRootFile{
    public:
        explicit PnpiRootFile(std::string& map_file);
        ~PnpiRootFile();
        void Init();

        TreeStructure* getTreeStructure();
        std::map<int,WaveFormParamStruct>* getSelfTreeStructure();
        std::set<int>* GetChannels();
        int GetEntries() const;
        uint GetSelfEntry(const uint ch);
        void GetNextEntry(int i);
        void GetNextEntry(uint i, uint ch);

    private:
        const std::string fileName_;
        TTree * allEvents_;
        TFile * fileInput_;
        std::set<int> availableChannels_;
        TreeStructure entry_;

        std::map<int, TTree *> allEventsSelfTrig_;
        std::map<int ,WaveFormParamStruct> sipmDataTrig_;
        int entries_;
    };
}

#endif //CUBESLY_PNPI_STRUCTURE_H
