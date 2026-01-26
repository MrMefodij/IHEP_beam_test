//
// Created by amefodev on 31.10.2024.
//

#include "FileHandler.h"
#include "defines.h"
#include <iostream>
#include <dirent.h>

std::vector<std::string> GetDataFiles(const std::string& stringBuf, const std::string& extension, const std::string& filter){
    std::vector<std::string> vFileNames;
    DIR *dir = opendir(stringBuf.c_str());
    if (dir == nullptr) {
        std::cerr << "Error: Cannot open directory " << stringBuf << std::endl;
    }
    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string fileName = entry->d_name;
        if (fileName.find(extension) != std::string::npos) {
            vFileNames.push_back(stringBuf + fileName);
        }
    }
    closedir(dir);
    auto it = std::remove_if(vFileNames.begin(), vFileNames.end(), [&](auto &item) {
        if (item.find(filter) != std::string::npos){
            std::cout << item << std::endl;
            return false;
        }
        return true;
    });
    vFileNames.erase(it, vFileNames.end());

    return vFileNames;
}

FileHandler::FileHandler(const std::string& inFileName, const std::string& extension, bool createOutput)
: inFileName_(inFileName), createOutput_(createOutput){
    if(!inFileName_.empty()) {
        size_t pos = inFileName_.rfind('.');
        if (pos != std::string::npos) {
            outFileName_ = inFileName_.substr(0, pos) + "_" + extension + ".root";
        }
        pos = inFileName_.rfind('/');
        if (pos != std::string::npos) {
            OutputDirName_ = inFileName_.substr(0, pos+1);
        }

        if (createOutput_) {
            wFile_ = std::make_shared<TFile>(outFileName_.c_str(), "recreate");
            if (!wFile_ || wFile_->IsZombie()) {
                std::cerr << "Failed to create output file: " << outFileName_ << std::endl;
                wFile_.reset();
            } else {
                std::cout << AS_KV(outFileName_) << " created." << "\n";
            }
        }
    }
}

std::shared_ptr<TFile> FileHandler::GetOutputFilePointer() {
    return wFile_;
}

std::string FileHandler::GetOutputDirName(){
    return OutputDirName_;
}

HitFileOutput::HitFileOutput(const std::string& inFileName,const std::string& extension) : FileHandler(inFileName, extension){
    EventsHits_ = std::make_shared<TTree>("EventsHits", "The ROOT tree of events");
    if (EventsHits_) {
        EventsHits_->Branch("Events", "PrototypeHit", &eventStructure_);
    } else {
        std::cerr << "Failed to create TTree EventsHits." << std::endl;
    }
}

PrototypeHit* HitFileOutput::GetEventPointer(){
    return &eventStructure_;
}

void HitFileOutput::FillTree(PrototypeHit& eventStructure){
    eventStructure_ = std::move(eventStructure);
    EventsHits_->Fill();
}

void HitFileOutput::Write(){
    wFile_->cd();
    EventsHits_->Write("", TObject::kOverwrite);
}


HitFileInput::HitFileInput(const std::string& inFileName,const std::string& extension, bool createOutput) : FileHandler(inFileName, extension, createOutput){
    if (!inFileName.empty()) {
        rFile_ = std::make_shared<TFile>(inFileName.c_str());
        if (!rFile_->IsOpen()) {
            std::cerr << "Can not open file " << inFileName << std::endl;
            return;
        }

        EventsHits_ = rFile_->Get<TTree>("EventsHits");

        if (!EventsHits_) {
            std::cerr << "Failed to retrieve TTree 'EventsHits' from file." << std::endl;
            return;
        }
        entries_ = EventsHits_->GetEntries();
        if (EventsHits_->GetBranch("Events")) {
            EventsHits_->SetBranchAddress("Events", &eventStructure_);
        }
        printf("File %s is opened.\n", inFileName.c_str());
    }
}

PrototypeHit* HitFileInput::getHitPointer(){
    return eventStructure_;
}

unsigned int HitFileInput::GetEntries() const{
    return entries_;
}

void HitFileInput::GetNextEntry(unsigned int i){
    EventsHits_->GetEntry(i);
}

void HitFileInput::Close() {
    if (rFile_) {
        rFile_->Close();
        rFile_.reset();
    }
}
