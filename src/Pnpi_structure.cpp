//
// Created by amefodev on 30.10.2024.
//

#include "Pnpi_structure.h"
#include <sstream>
#include <iomanip>
#include <TKey.h>

namespace PNPI{
    PrototypeHit GetHitStructure (const TreeStructure* entry, const std::set<int>* availableChannels,const Calibration* calibration) {
        PrototypeHit result;

        if (!entry || !availableChannels || !calibration) {
            std::cerr << "ERROR: One of the pointers is NULL in GetHitStructure!" << std::endl;
            return result;
        }

        const std::map<int, CalibResults> *calibrationResultsI = calibration->GetCalibResultsI();
        const std::map<int, CalibResults> *calibrationResultsA = calibration->GetCalibResultsA();
        const std::map<int, bool> *channelsFlag = calibration->GetChannelFlag();
        result._quality = (int) entry->_trackParam.Q;
        result._x = entry->_trackParam.AX * CUBES_POSITION + entry->_trackParam.BX;
        result._y = entry->_trackParam.AY * CUBES_POSITION +
                    entry->_trackParam.BY; //cube position 1450 +/- 1 cm (1485) was 1450 become 1485

        result._Ax = entry->_trackParam.AX;
        result._Bx = entry->_trackParam.BX;
        result._Ay = entry->_trackParam.AY;
        result._By = entry->_trackParam.BY;

        for (const auto ch: *availableChannels) {
            if (entry->_sipmData.count(ch) == 0) continue;
            if (ch == -1) continue;

            if (entry->_sipmData.at(ch).A > 0) {
                ChannelData hit;

                if (channelsFlag != nullptr && channelsFlag->count(ch)) {
                    hit._channels_flag = channelsFlag->at(ch);
                } else {
                    hit._channels_flag = false;
                }

                hit._channel_id = ch;
                hit._quality = entry->_sipmData.at(ch).Q;
                hit._amplitude = entry->_sipmData.at(ch).A;
                hit._time = entry->_sipmData.at(ch).T;
                hit._chargeADC = entry->_sipmData.at(ch).I;

                if (calibrationResultsI != nullptr){
                    hit._chargePE_Integral = (entry->_sipmData.at(ch).I - calibrationResultsI->at(ch)._peak0) / calibrationResultsI->at(ch)._gain + 1;
//                    hit._chargePE_Integral_XTalk = (entry->_sipmData.at(ch).I - calibrationResultsI->at(ch)._peak0) / (calibrationResultsI->at(ch)._gain * calibrationResultsI->at(ch)._xTalk) + 1;
                } else {
                    std::cerr << "Сhannel ID: " << ch << " Нет калибровки I " << ch << std::endl;
                    hit._chargePE_Integral = 0;
                }
                if (calibrationResultsA != nullptr){
                    hit._chargePE_Amplitude = (entry->_sipmData.at(ch).A - calibrationResultsA->at(ch)._peak0) / calibrationResultsA->at(ch)._gain + 1;
//                    hit._chargePE_Amplitude_XTalk = (entry->_sipmData.at(ch).A - calibrationResultsA->at(ch)._peak0) / (calibrationResultsA->at(ch)._gain * calibrationResultsA->at(ch)._xTalk) + 1;
                } else {
                    std::cerr << "Сhannel ID: " << ch << " Нет калибровки A " << ch << std::endl;
                    hit._chargePE_Integral = 0;
                }

                result._channelsData.push_back(hit);
            }
        }
        return result;
    }

    PnpiRootFile::PnpiRootFile(std::string& file) : fileName_(file){}

    PnpiRootFile::~PnpiRootFile() {
        delete allEvents_;
        fileInput_->Close();
    }

    void PnpiRootFile::Init(){
        fileInput_ = new TFile(fileName_.c_str());
        if (!fileInput_->IsOpen()) {
            std::cerr << "Can not open file " << fileName_ << std::endl;
            return;
        }
        TIter nextKey(fileInput_->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)nextKey())) {
            TString keyName = key->GetName();
            if (keyName.BeginsWith("Event")) {
                allEvents_ = (TTree *) fileInput_->Get(keyName);
                entries_ = allEvents_->GetEntries();
                if (allEvents_->GetBranch("Track")) {
                    allEvents_->SetBranchAddress("Track", &entry_._trackParam);
                }

                std::ostringstream sChNum;
                std::string sCh;
                for (int ih = 0; ih <= CHANNELS_NUMBER; ++ih) {
                    sChNum.str("");
                    sChNum << std::setw(2) << std::setfill('0') << ih;
                    sCh = "SiPM_55_" + sChNum.str();
                    if (allEvents_->GetBranch((sCh).c_str())) {
                        entry_._sipmData[ih] = {};
                        availableChannels_.insert(ih);
                        allEvents_->SetBranchAddress((sCh).c_str(), &entry_._sipmData[ih]);
                    }
                }

                for (int ih = 0; ih <= CHANNELS_NUMBER; ++ih) {
                    sChNum.str("");
                    sChNum << std::setw(2) << std::setfill('0') << ih;
                    sCh = "SiPM_55_" + sChNum.str();
                    allEventsBeamTrig_[ih] = (TTree*) (fileInput_->Get(keyName));
                    if (allEventsBeamTrig_[ih]->GetBranch((sCh).c_str())) {
                        beamDataTrig_[ih] = {};
                        allEventsBeamTrig_[ih]->SetBranchAddress((sCh).c_str(), &beamDataTrig_[ih]);
                    }
                }
            }
            else {
                std::string s = "Self_55_";
                std::string rootTree = static_cast<std::string>(keyName);
                if (keyName.BeginsWith(s.c_str())) {
                     auto ih = std::stoul(rootTree.substr(rootTree.find(s) + s.length(), 2));
                    allEventsSelfTrig_[ih] = (TTree*) (fileInput_->Get(keyName));
                    if (allEventsSelfTrig_[ih]->GetBranch("Self")) {
                        sipmDataTrig_[ih] = {};
                        allEventsSelfTrig_[ih]->SetBranchAddress("Self", &sipmDataTrig_[ih]);
                    }
                }
            }
        }
    }

    TreeStructure* PnpiRootFile::getTreeStructure(){
        return &entry_;
    }

    std::map<int,WaveFormParamStruct>* PnpiRootFile::getSelfTreeStructure(){
        return &sipmDataTrig_;
    }

    std::map<int,WaveFormParamStruct>* PnpiRootFile::getBeamTreeStructure(){
        return &beamDataTrig_;
    }

    std::set<int>* PnpiRootFile::GetChannels(){
        return &availableChannels_;
    };

    int PnpiRootFile::GetEntries() const{
        return entries_;
    }

    uint PnpiRootFile::GetSelfEntry(const uint ch){
        auto iter = allEventsSelfTrig_.find(ch);
        if(iter !=  allEventsSelfTrig_.end())
            return allEventsSelfTrig_[ch]->GetEntries();
        else{
            std::cerr << "Channel \'" << ch << "\' isn't appropriate" << "\'\n";
        }
        return 0.0;
    }

    uint PnpiRootFile::GetBeamEntry(const uint ch){
        auto iter = allEventsBeamTrig_.find(ch);
        if(iter !=  allEventsBeamTrig_.end())
            return allEventsBeamTrig_[ch]->GetEntries();
        else{
            std::cerr << "Channel \'" << ch << "\' isn't appropriate" << "\'\n";
        }
        return 0.0;
    }

    void PnpiRootFile::GetNextEntry(int i){
        if (i <= entries_) {
            allEvents_->GetEntry(i);
        } else {
            std::cerr << "Entry \'" << i << "\' is more than number of entries \'" << entries_ << "\'\n";
        }
    }
    void PnpiRootFile::GetNextEntry(uint i, uint ch){
        allEventsSelfTrig_.at(ch)->GetEntry(i);
    }
    void PnpiRootFile::GetNextBeamEntry(uint i, uint ch){
        allEventsBeamTrig_.at(ch)->GetEntry(i);
    }
}