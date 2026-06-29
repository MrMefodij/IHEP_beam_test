//
// Created by Angelina Chvirova on 22.06.2026.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <TF1.h>
#include "TFile.h"
#include "TH1F.h"
#include "connection_map.h"
#include "FileHandler.h"

int main(int argc, char **argv) {
    std::string filename = argv[1];
    Connection_Map *connectionMap = new Connection_Map();
    connectionMap->Init();
    const std::map<int, double> *fibers_normalization = connectionMap->GetFibersNormalization();

    if (!fibers_normalization) {
        return 1;
    }

    std::map<int, TH1F *> hists;
    TH1F *h_all_SFGD = new TH1F("h_all_SFGD", "LY All MPPC SFGD Channels; LY [p.e.]; N", FIBER_BINNING);
    TH1F *h_all_BabyM = new TH1F("h_all_BabyM", "LY All MPPC Baby Mind Channels; LY [p.e.]; N", FIBER_BINNING);

    size_t pos = filename.rfind('/');
    std::string OutputDirName = (pos != std::string::npos) ? filename.substr(0, pos + 1) : "";

    std::shared_ptr<TFile> out = std::make_shared<TFile>((OutputDirName + "LY_Fibers.root").c_str(), "recreate");
    std::vector<std::string> dataFiles = GetDataFiles(OutputDirName, ".root", "_HitStructure.");


    for (const auto &file : dataFiles) {
        size_t lastSlash = file.rfind('/');
        std::string fileNameOnly = (lastSlash == std::string::npos) ? file : file.substr(lastSlash + 1);

        if (fileNameOnly.empty() || fileNameOnly[0] == '.' || fileNameOnly.find("._") != std::string::npos) {
            continue;
        }


        try {
            HitFileInput files(file, "cubes_position", false);
            int entries = files.GetEntries();

            for (int i = 0; i < entries; ++i) {
                files.GetNextEntry(i);
                auto *hit = files.getHitPointer();

                if (!hit) continue;


                for (const auto &ch : hit->_channelsData) {
                    int id = ch._channel_id;
                    if (id < 0) continue;

                    if (ch._chargePE_Amplitude > 1 && ch._quality < 100) {

                        auto itNorm = fibers_normalization->find(id);
                        if (itNorm == fibers_normalization->end() || itNorm->second == 0) {
                            continue;
                        }

                        double ly = ch._chargePE_Amplitude / itNorm->second;


                        if (hists.find(id) == hists.end()) {
                            hists[id] = new TH1F(Form("h_ch_%d", id), Form("LY Channel %d; LY [p.e.]; N", id), FIBER_BINNING);
                            hists[id]->SetDirectory(0);
                        }

                        hists[id]->Fill(ly);

                        if (id < 50) {
                            h_all_SFGD->Fill(ly);
                        } else {
                            h_all_BabyM->Fill(ly);
                        }
                    }
                }
            }
        } catch (const std::exception& e) {
            continue;
        }
    }


    out->cd();
    for (auto &[id, hist] : hists) {
        if (!hist || hist->GetEntries() <= 50) continue;

        if (id < 50) {
            TF1 *fit = new TF1(Form("gaus_ch_%d", id), "gaus", 15, 25);
            hist->Fit(fit, "qr+");
            hist->Write();
            delete fit;
        } else if (id > 50){
            TF1 *fit = new TF1(Form("gaus_ch_%d", id), "gaus", 25, 55);
            hist->Fit(fit, "qr+");
            hist->Write();
            delete fit;
        }
    }

    if (h_all_SFGD->GetEntries() > 50) {
        TF1 *fitSFGD = new TF1("fit_all_SFGD", "gaus", 15, 25);
        h_all_SFGD->Fit(fitSFGD, "qr+");
        h_all_SFGD->Write();
        delete fitSFGD;
    }

    if (h_all_BabyM->GetEntries() > 50) {
        TF1 *fitBabyM = new TF1("fit_all_BabyM", "gaus", 25, 55);
        h_all_BabyM->Fit(fitBabyM, "qr+");
        h_all_BabyM->Write();
        delete fitBabyM;
    }

    out->Close();

    delete connectionMap;
    for (auto &[id, hist] : hists) delete hist;
    delete h_all_SFGD;
    delete h_all_BabyM;
    return 0;
}