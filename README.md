LFVStatistics_v3
================

The macro Sens_discovery.C computes the 3sigma sensitivity for discovery using "data" and "signal" root files
containing the EM and ME collinear mass histograms.
****************************
Setup before running:
  in run/LoadAll.C:
    change the following line to include your own path (where this directory sits in your file system):
    
    gSystem->AddIncludePath("-I/home/avitald/LFV/LFVStatistics_v3");
    
  in Root/Constants.cxx:
    change the definition of the histogam names to fit the naming in your root files:
    
     HISTO_NAME_ME = "nom/ME_McollHiggs_Unblind";
     HISTO_NAME_EM = "nom/EM_McollHiggs_Unblind";

****************************
Running the macro:
  
  cd run/
  .x LoadAll.C
  Sens_discovery(data_root_file_path,signal_root_file_path)
  
****************************

Output:
1) A "brazil plot" of the sensitivity histogram.
2) the median sensitivity value and 1/2sigma band values are printed to the screen.
3) the values are also written to a file called SensFile.txt
