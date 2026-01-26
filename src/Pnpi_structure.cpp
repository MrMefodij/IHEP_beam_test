//
// Created by amefodev on 30.10.2024.
//

#include "Pnpi_structure.h"
#include <sstream>
#include <iomanip>
#include <TKey.h>

namespace PNPI{
    PrototypeHit GetHitStructure (const TreeStructure* entry, const std::set<unsigned int>* availableChannels,const Calibration* calibration){
        PrototypeHit result;
        const std::map<unsigned int,CalibResults>* calibrationResultsI = calibration->GetCalibResultsI();
        const std::map<unsigned int,CalibResults>* calibrationResultsA = calibration->GetCalibResultsA();
        const std::map<unsigned int, bool>* channelsFlag = calibration->GetChannelFlag();
        result._quality = (int)entry->_trackParam.Q;
        result._x = entry->_trackParam.AX * CUBES_POSITION + entry->_trackParam.BX;
        result._y = entry->_trackParam.AY * CUBES_POSITION + entry->_trackParam.BY; //cube position 1450 +/- 1 cm (1485) was 1450 become 1485

        result._Ax = entry->_trackParam.AX;
        result._Bx = entry->_trackParam.BX;
        result._Ay = entry->_trackParam.AY;
        result._By = entry->_trackParam.BY;

        for (const auto ch : *availableChannels) {
            if (entry->_sipmData.at(ch).A > 0){
                ChannelData hit;
                if (channelsFlag != nullptr) {
                    hit._channels_flag = channelsFlag->at(ch);
                }
                hit._channel_id = ch ;
                hit._quality = entry->_sipmData.at(ch).Q;
                hit._amplitude = entry->_sipmData.at(ch).A;
                hit._time = entry->_sipmData.at(ch).T;
                hit._chargeADC = entry->_sipmData.at(ch).I;
                hit._reconstructed_amplitude = entry->_sipmData.at(ch).AF;
                if (calibrationResultsI != nullptr){
                    hit._chargePE_Integral = (entry->_sipmData.at(ch).I - calibrationResultsI->at(ch)._peak0) / calibrationResultsI->at(ch)._gain + 1;
                    hit._chargePE_Integral_XTalk_with_background = (entry->_sipmData.at(ch).I - calibrationResultsI->at(ch)._peak0) / (calibrationResultsI->at(ch)._gain * calibrationResultsI->at(ch)._xTalk_with_background) + 1;
                    hit._chargePE_Integral_XTalk_without_background = (entry->_sipmData.at(ch).I - calibrationResultsI->at(ch)._peak0) / (calibrationResultsI->at(ch)._gain * calibrationResultsI->at(ch)._xTalk_without_background) + 1;
                }
                if (calibrationResultsA != nullptr){
                    hit._chargePE_Amplitude = (entry->_sipmData.at(ch).A - calibrationResultsA->at(ch)._peak0) / calibrationResultsA->at(ch)._gain + 1;
                    hit._chargePE_Amplitude_XTalk_with_background = (entry->_sipmData.at(ch).A - calibrationResultsA->at(ch)._peak0) / (calibrationResultsA->at(ch)._gain * calibrationResultsA->at(ch)._xTalk_with_background) + 1;
                    hit._chargePE_Amplitude_XTalk_without_background = (entry->_sipmData.at(ch).A - calibrationResultsA->at(ch)._peak0) / (calibrationResultsA->at(ch)._gain * calibrationResultsA->at(ch)._xTalk_without_background) + 1;
                    hit._reconstructed_chargePE_Amplitude = (entry->_sipmData.at(ch).AF - calibrationResultsA->at(ch)._peak0) / (calibrationResultsA->at(ch)._gain * calibrationResultsA->at(ch)._xTalk_with_background) + 1;
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
                for (unsigned int ih = 0; ih <= CHANNELS_NUMBER; ++ih) {
                    sChNum.str("");
                    sChNum << std::setw(2) << std::setfill('0') << ih;
                    sCh = "SiPM_55_" + sChNum.str();
                    if (allEvents_->GetBranch((sCh).c_str())) {
                        entry_._sipmData[ih] = {};
                        availableChannels_.insert(ih);
                        allEvents_->SetBranchAddress((sCh).c_str(), &entry_._sipmData[ih]);
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

    std::map<unsigned int,WaveFormParamStruct>* PnpiRootFile::getSelfTreeStructure(){
        return &sipmDataTrig_;
    }

    std::set<unsigned int>* PnpiRootFile::GetChannels(){
        return &availableChannels_;
    };

    unsigned int PnpiRootFile::GetEntries() const{
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

    void PnpiRootFile::GetNextEntry(unsigned int i){
        if (i <= entries_) {
            allEvents_->GetEntry(i);
        } else {
            std::cerr << "Entry \'" << i << "\' is more than number of entries \'" << entries_ << "\'\n";
        }
    }
    void PnpiRootFile::GetNextEntry(uint i, uint ch){
        allEventsSelfTrig_.at(ch)->GetEntry(i);
    }
}