//
// Created by Angelina Chvirova on 17.08.2026.
//

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <iostream>

int main() {
//    const int nPoints = 16; // 150 - 168
//    double x[nPoints] = {150.0, 151.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 168.0};
//    double y[nPoints] = {19.2, 19.19,19.14, 19.12, 19.09, 19.0,  18.94, 18.9, 18.85, 18.78, 18.75, 18.75, 18.73, 18.71, 18.7 , 18.63};

    const int nPoints = 14; // 96 - 110
    double x[nPoints] = {96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0};
    double y[nPoints] = {24.57, 24.49,24.46, 24.45, 24.42, 24.41,  24.46, 18.9, 18.85, 18.78, 18.75, 18.75, 18.73, 18.71, 18.7 , 18.63};


    TGraph *gr = new TGraph(nPoints, x, y);
    gr->SetName("LY");
    gr->SetTitle("LY run 150 – 168; Run number; LY [p.e.]");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerSize(1.2);

//    TF1 *fitFunc = new TF1("fitFunc", "pol1", x[0], x[nPoints - 1]); // "pol1" -> p0 + p1*x
//    fitFunc->SetLineColor(kRed);
//    fitFunc->SetLineWidth(2);
//    gr->Fit(fitFunc, "R");

    TFile *outFile = new TFile("LY_graphs.root", "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file for writing!" << std::endl;
        return 1;
    }

    // Записываем граф и функцию в файл
    gr->Write("LY_");
//    fitFunc->Write("Fit_LY");

    outFile->Close();
    delete outFile;

//    TCanvas *c1 = new TCanvas("c1", "Linear Fit", 800, 600);
//    gr->Draw("AP");
//    c1->Modified();
//    c1->Update();
    return 0; // Обязательно для main()
}