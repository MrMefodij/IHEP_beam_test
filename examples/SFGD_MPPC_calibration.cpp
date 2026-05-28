//
// Created by Maria on 07.11.2024.
//
#include <string>
#include <TH1I.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <fstream>
#include <TGraphErrors.h>
#include <TDirectory.h>
#include "connection_map.h"
#include "Pnpi_structure.h"
#include "FileHandler.h"

void Calib(const std::string& filename, PNPI::PnpiRootFile* fileInput,  FileHandler* fileOut, TDirectory* directory, const char charge) {
    {
        TH1F gain = TH1F("gain","gain",100,0,20);
        TH1F x_talk = TH1F("x_talk","x_talk",100,0, 30);

        std::string fileOutput;
        size_t pos = filename.rfind('.');
        if (pos != std::string::npos) {
            if (charge == 'A') {
                fileOutput = filename.substr(0, pos)  + "_Calibration_A.txt";
            } else {
                fileOutput = filename.substr(0, pos)  + "_Calibration_I.txt";
            }
        }
        std::ofstream outFile(fileOutput);

        if (!outFile) {
            std::cerr << "Error while openning file: "<< fileOutput << std::endl;
            return;
        }
        outFile <<"ch gain gain_error peak_1 peak_1_error mean_w/_background x-talk_w/_background mean_w/o_background x-talk_w/o_background\n";

        std::map<int , PNPI::WaveFormParamStruct>* hits = fileInput->getSelfTreeStructure();
        std::map<int , PNPI::WaveFormParamStruct>* hits_beam = fileInput->getBeamTreeStructure();
        std::set<int>* availableChannels = fileInput->GetChannels();

        TSpectrum *fTSpectrum = new TSpectrum(60);
        TF1 *f1 = new TF1("Linear fit", "[0]*x+[1]");
        f1->SetLineColor(kRed);
        f1->SetLineStyle(2);
        f1->SetLineWidth(4);
        fileOut->Write();
        directory->cd();
        for (auto  ch : *availableChannels) {
            std::string sCh = "Channel_" + std::to_string(ch);
            auto entry = fileInput->GetSelfEntry(ch);
            if (fileInput->GetSelfEntry(ch) == 0){ continue;}
            std::cout << "Getting channel " << ch << " with " << entry << " number of enters." << std::endl;
            TH1I * hFEBCH;
            TH1I * hCH_xtalk;
            if (charge == 'A') {
                hFEBCH = new TH1I("hFEBCH","hFEBCH", 70,0,70);
                hCH_xtalk = new TH1I("hCH_xtalk","hCH_xtalk", 700,0,700);
            } else {
                hFEBCH = new TH1I("hFEBCH","hFEBCH", 300,0,600);
                hCH_xtalk = new TH1I("hCH_xtalk","hCH_xtalk", 300,0,600);
            }

            TGraphErrors *gr = new TGraphErrors;
            for (int i = 0; i < entry; ++i) {
                fileInput->GetNextEntry(i, ch);
                if  (hits->at(ch).Q < 1000) {
                    if (hits->at(ch).A > 5) {
                        if (charge == 'A') {
                            hFEBCH->Fill(hits->at(ch).A);
                        } else {
                            hFEBCH->Fill(hits->at(ch).I);
                        }
                        if (hits->at(ch).Q <= 2) {
                            if (charge == 'A') {
                                hCH_xtalk->Fill(hits->at(ch).A);
                            } else {
                                hCH_xtalk->Fill(hits->at(ch).I);
                            }
                        }
                    }
                }
            }

            auto beam_entry = fileInput->GetBeamEntry(ch);
            for (int i = 0; i < beam_entry; ++i) {
                fileInput->GetNextBeamEntry(i, ch);
                if  (hits_beam->at(ch).Q < 1000) {
                    if (hits_beam->at(ch).A > 5) {
                        if (charge == 'A') {
                            hFEBCH->Fill(hits_beam->at(ch).A);
                        } else {
                            hFEBCH->Fill(hits_beam->at(ch).I);
                        }
                    }
                }
            }

            double mean_w_background = hCH_xtalk->GetMean();
            hFEBCH->SetTitle((sCh + "_Calibration" ).c_str());
            hFEBCH->SetName((sCh + "_Calibration" ).c_str());

            hCH_xtalk->SetTitle((sCh + "_X-Talk" ).c_str());
            hCH_xtalk->SetName((sCh + "_X-Talk" ).c_str());

            if (charge == 'A') {
                hFEBCH->GetXaxis()->SetTitle("Amplitude [ADC]");
            } else {
                hFEBCH->GetXaxis()->SetTitle("Integral [ADC]");
            }
            hFEBCH->GetYaxis()->SetTitle("N");
            hCH_xtalk->GetYaxis()->SetTitle("N");
            hCH_xtalk->Write();

            TH1* histBackground;
            if (charge == 'A') {
                histBackground = fTSpectrum->Background(hFEBCH, 10, "");
            } else {
                histBackground = fTSpectrum->Background(hFEBCH, 20, "");
            }
            // hFEBCH->Add(histBackground, -1);

            int nfound;
            if (charge == 'A') {
                nfound = fTSpectrum->Search(hFEBCH,2,"",0.001);
            } else {
                nfound = fTSpectrum->Search(hFEBCH,4,"",0.002);
            }

            if(nfound > 0) {
                double *xpeaks = fTSpectrum->GetPositionX();
                std::sort(xpeaks, xpeaks + nfound);
                double pedestal {0.0};
                double pedestalError {0.0};
                uint peakNumber = 0;
                for (uint i = 0; i < nfound; ++i) {
                    double first_peak;
                    double cut = 40;
                    if (charge == 'A') {
                        first_peak = 6.5;
                    } else {
                        first_peak = 40;
                    }
                    double last_peak = xpeaks[0];
                    if ( (i == 0 && xpeaks[0] > first_peak) || ( i != 0  && xpeaks[i] > xpeaks[i - 1] && xpeaks[i] > last_peak + 5 && xpeaks[i] < cut) )  {
                        TF1* fit_1;
                        if (charge == 'A') {
                            fit_1 = new TF1("fit_1", "gaus", xpeaks[i] - 3, xpeaks[i] + 3);
                        } else {
                            fit_1 = new TF1("fit_1", "gaus", xpeaks[i] - 10, xpeaks[i] + 10);
                        }
                        last_peak = xpeaks[i];

                        hFEBCH->Fit("fit_1", "qr+");
                        if (peakNumber == 0) {
                            pedestal = fit_1->GetParameter(1);
                            pedestalError = fit_1->GetParameter(2);
                        }
                        gr->AddPoint(peakNumber + 1, fit_1->GetParameter(1));
                        gr->SetPointError(peakNumber, 0, fit_1->GetParameter(2));
                        delete fit_1;
                        peakNumber++;
                    }
                }
                gr->Fit(f1, "qp");
                // pedestal = pedestal - f1->GetParameter(0);
                double x_talk_results = (mean_w_background - (pedestal -  f1->GetParameter(0)))/ f1->GetParameter(0);
                if (x_talk_results < 1) {
                    x_talk_results = 1;
                }
                outFile << ch << " " << f1->GetParameter(0) << " "<<  f1->GetParError(0) <<" "<< pedestal<< " "<< pedestalError << " " <<
                    mean_w_background << " " << x_talk_results << "\n";
                gain.Fill(f1->GetParameter(0));
                x_talk.Fill((x_talk_results) * 100 - 100);

                gr->SetTitle(("Fit: Channel_"+ std::to_string(ch)).c_str());
                gr->SetName(("Fit: Channel_"+ std::to_string(ch)).c_str());
                gr->GetXaxis()->SetTitle("Peak number");
                gr->GetYaxis()->SetTitle("ADC channels");
                gr->Write();
            }
            histBackground->Delete();
            hFEBCH->Write();
            hFEBCH->Delete();
            gr->Delete();
        }
        gain.GetXaxis()->SetTitle("ADC");
        gain.GetYaxis()->SetTitle("N");
        gain.Write();
        x_talk.GetXaxis()->SetTitle("%");
        x_talk.GetYaxis()->SetTitle("N");
        x_talk.Write();

        fTSpectrum->Delete();
        f1->Delete();
        outFile.close();
    }
}

int main( int argc, char **argv ) {
    std::string filename = argv[1];

    Connection_Map connectionMap;
    connectionMap.Init();

    PNPI::PnpiRootFile fileInput(filename);
    fileInput.Init();

    FileHandler fileOut(filename, "Calibration");
    std::shared_ptr<TFile> rFile = fileOut.GetOutputFilePointer();
    TDirectory* Ampl = rFile->mkdir("Amplitude");
    // TDirectory* Int = rFile->mkdir("Integral");
    Calib(filename, &fileInput,  &fileOut, Ampl, 'A');
    // Calib(filename, &fileInput,  &fileOut, Int, 'I');
    return 0;
}