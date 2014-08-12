#include "TROOT.h"
#include "TSystem.h"

void LoadAll()
{

    gROOT->ProcessLine("#include <vector>");

    // Set the include path for all packages
    gSystem->AddIncludePath("-I/home/avitald/LFV/LFVStatistics_v3");
    cout << gSystem->GetIncludePath() << endl;

    // Load package libraries
    gROOT->ProcessLine(".L ../cmt/LFVStatistics_v3.C+");

    //Load additional macros
    gROOT->ProcessLine(".L ../macros/Start.C+");
    gROOT->ProcessLine(".L ../macros/Sens_discovery.C+");
    gROOT->ProcessLine(".L ../macros/pdf_limit.C+");
    gROOT->ProcessLine(".L ../macros/pdf_discovery.C+");


}


