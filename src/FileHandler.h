//
// Created by amefodev on 31.10.2024.
//

#ifndef CUBESLY_FILEHANDLER_H
#define CUBESLY_FILEHANDLER_H
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "PrototypeHit.h"

ClassImp(ChannelData)
ClassImp(PrototypeHit)

std::vector<std::string> GetDataFiles(const std::string& stringBuf, const std::string& extension, const std::string& filter);

class FileHandler {
public:
    explicit FileHandler(const std::string& inFileName, const std::string& extension, bool createOutput = true);
    std::shared_ptr<TFile> GetOutputFilePointer();
    std::string GetOutputDirName();
    void Write(){
        wFile_->cd();
    }
private:
    const std::string inFileName_;
    const bool createOutput_;
    std::string outFileName_;
    std::string OutputDirName_;
protected:
    std::shared_ptr<TFile> wFile_ = nullptr;
};

class HitFileOutput : public FileHandler{
public:
    explicit HitFileOutput(const std::string& inFileName,const std::string& extension);
    PrototypeHit* GetEventPointer();
    void FillTree(PrototypeHit& eventStructure);
    void Write();
private:
    std::shared_ptr<TTree> EventsHits_ = nullptr;
    PrototypeHit eventStructure_;
};

class HitFileInput : public FileHandler{
public:
    explicit HitFileInput(const std::string& inFileName,const std::string& extension, bool createOutput = true);
    PrototypeHit* getHitPointer();
    unsigned int GetEntries() const;
    void GetNextEntry(unsigned int i);
    void Close();

private:
    std::shared_ptr<TFile> rFile_ = nullptr;
    unsigned int entries_;
    TTree* EventsHits_ = nullptr;
    PrototypeHit* eventStructure_;
};



#endif //CUBESLY_FILEHANDLER_H
