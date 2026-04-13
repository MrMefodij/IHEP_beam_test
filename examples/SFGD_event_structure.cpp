//
// Created by amefodev on 30.10.2024.
//

#include <string>
#include "defines.h"
#include "Calibration.h"
#include "Pnpi_structure.h"
#include "FileHandler.h"
#include "PrototypeHit.h"

int main( int argc, char **argv ) {
    std::string filename = argv[1];

    PNPI::PnpiRootFile fileInput(filename);
    fileInput.Init();
    Calibration calibration(filename);

    HitFileOutput fileOut(filename, "HitStructure");
    std::shared_ptr<TFile> rFile = fileOut.GetOutputFilePointer();

    PNPI::TreeStructure* hits = fileInput.getTreeStructure();

    if (!hits) {
        std::cerr << "Error: TreeStructure is null!" << std::endl;
        return 1;
    }

    std::set<int>* availableChannels = fileInput.GetChannels();

    for (int i = 0; i < fileInput.GetEntries(); ++i) {
        fileInput.GetNextEntry(i);
        PrototypeHit events = PNPI::GetHitStructure(hits,availableChannels, &calibration);
        if (events._quality>=0 && events._quality<=20) {
            fileOut.FillTree(events);
        }
    }
    fileOut.Write();
    return 0;
}