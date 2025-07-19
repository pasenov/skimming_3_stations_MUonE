// Converting MC simulated files (only simulation and digitization) to the same format as data
// Exception: for simulated data we added MC branches

#include <array>
#include <cmath>      // for std::abs and std::fabs
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "FairMCEventHeader.h"
#include "MUonETrack.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TFolder.h"
#include "TLegend.h"
#include "TList.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TVector3.h"


const int SAVE = 1;
const int NMODULES = 12;

void runConvertMC(TString filename)
{
    TString filepath = "/eos/user/p/pasenov/www/pasenov/MUonE_data/Condor_results/3stations_2025/MCminBias_3stations/";

    // TString filepath = "/eos/user/p/pasenov/www/pasenov/MUonE_data/Condor_results/3stations_2025/MCsignalJob_3stations/";
    TString fullpath = filepath + filename;

    TFile *infile = new TFile(fullpath);
    TTree *cbmsim_input = (TTree *)infile->Get("cbmsim");

    TClonesArray *TrackerStubs = nullptr;
    TClonesArray *MCTrack = nullptr;
    FairMCEventHeader *MCEventHeader = nullptr;
    MUonERecoOutputFull *ReconstructionOutput = nullptr;
    MCEventTruth *MCEventTruth = nullptr;

    cbmsim_input->SetBranchAddress("TrackerStubs", &TrackerStubs);
    cbmsim_input->SetBranchAddress("MCTrack", &MCTrack);
    cbmsim_input->SetBranchAddress("MCEventHeader.", &MCEventHeader);
    cbmsim_input->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);
    cbmsim_input->SetBranchAddress("MCEventTruth", &MCEventTruth);

    TString outputfilename = fullpath;
    outputfilename.Remove(outputfilename.Sizeof() - 1 - 5);
    outputfilename += "_convertedTFormat.root";
    TFile *outfile = new TFile(outputfilename, "RECREATE");
    TTree *output_cbmsim = new TTree("cbmsim", "cbmsim");

    // Copy the branches from the input tree to the output tree
    output_cbmsim->Branch("MCTrack", &MCTrack);
    output_cbmsim->Branch("TrackerStubs", &TrackerStubs);
    output_cbmsim->Branch("MCEventHeader", &MCEventHeader);
    output_cbmsim->Branch("ReconstructionOutput", &ReconstructionOutput);
    output_cbmsim->Branch("MCEventTruth", &MCEventTruth);

    // New branches
    std::vector<unsigned short> Bx;
    std::vector<unsigned short> Link;
    std::vector<float> LocalX;
    std::vector<float> LocalY;
    std::vector<float> Bend;
    std::vector<unsigned int> SuperID;
    std::vector<unsigned short> StationID;


    output_cbmsim->Branch("Bx", &Bx);
    output_cbmsim->Branch("Link", &Link);
    output_cbmsim->Branch("LocalX", &LocalX);
    output_cbmsim->Branch("LocalY", &LocalY);
    output_cbmsim->Branch("Bend", &Bend);
    output_cbmsim->Branch("SuperID", &SuperID);
    output_cbmsim->Branch("StationID", &StationID);

    TList *tl = new TList();
    tl->SetOwner();
    tl->AddLast(new TObjString("MCEventTruth"));
    tl->AddLast(new TObjString("ReconstructionOutput"));
    tl->AddLast(new TObjString("MCTrack"));
    tl->AddLast(new TObjString("MCEventHeader"));
    tl->AddLast(new TObjString("TrackerStubs"));
    tl->AddLast(new TObjString("Bx"));
    tl->AddLast(new TObjString("Link"));
    tl->AddLast(new TObjString("LocalX"));
    tl->AddLast(new TObjString("LocalY"));
    tl->AddLast(new TObjString("Bend"));
    tl->AddLast(new TObjString("SuperID"));
    tl->AddLast(new TObjString("StationID"));
    tl->Write("BranchList", 1);
    delete tl;

    // TH1F *h_moduleID = new TH1F("h_module", "ModuleID ", 8, -1, 7);

    unsigned int superid = 0;
    unsigned short bx = 0;
    int counter = 0;

    for (int i = 0; i < cbmsim_input->GetEntries(); i++)
    {
        cbmsim_input->GetEntry(i);

        if (i % 50000 == 0)
            std::cout << "\rprocessing event : " << i << std::flush << endl;

        bx = i % 3564;
        if (i % 3564 == 0)
            superid++;
        if (i == cbmsim_input->GetEntries() - 1)
        {
            std::cout << "Last superid: " << superid << std::endl;
        }

        if (TrackerStubs->GetEntries() == 0)
            continue;

        for (int j = 0; j < TrackerStubs->GetEntries(); j++)
        {
            const MUonETrackerStub *MCStub = static_cast<MUonETrackerStub *>(TrackerStubs->At(j));

            unsigned short link = MCStub->stationID() * 6 + MCStub->moduleID();
            // h_moduleID->Fill(MCStub->moduleID());
            unsigned short stationID = MCStub->stationID();

            float localx = MCStub->seedClusterCenterStrip();
            float bend = MCStub->bend();
            if (MCStub->moduleID() % 6 == 0 || MCStub->moduleID() % 6 == 4)
            {
                localx = 1015 - MCStub->seedClusterCenterStrip();
                bend = -bend;
            }
            float localy = MCStub->sensorHalf() == 0 ? 0.25 : 0.75;

            Bx.push_back(bx);
            Link.push_back(link);
            LocalX.push_back(localx);
            LocalY.push_back(localy);
            Bend.push_back(bend);
            SuperID.push_back(superid);
            StationID.push_back(stationID);
        }

        output_cbmsim->Fill();
        counter++;

        Bx.clear();
        Link.clear();
        LocalX.clear();
        LocalY.clear();
        Bend.clear();
        SuperID.clear();
        StationID.clear();
    }

    output_cbmsim->Write();

    std::cout << std::endl;
    std::cout << "nevents written = " << counter << std::endl;

    // Recreate the TFolder
    TFolder *originalFolder = (TFolder *)infile->Get("cbmroot");
    TFolder *newFolder = (TFolder *)originalFolder->Clone();
    outfile->Add(newFolder);

    // Recreate the TimeBasedBranchList
    TList *timeBasedBranchList = (TList *)infile->Get("TimeBasedBranchList");
    TList *newTimeBasedBranchList = (TList *)timeBasedBranchList->Clone();
    newTimeBasedBranchList->SetName("TimeBasedBranchList");
    outfile->Add(newTimeBasedBranchList);

    // Recreate the FileHeader
    FairFileHeader *fileHeader = (FairFileHeader *)infile->Get("FileHeader");
    FairFileHeader *newFileHeader = (FairFileHeader *)fileHeader->Clone();
    newFileHeader->SetName("FileHeader");
    outfile->Add(newFileHeader);

    outfile->cd();
    outfile->Write();


    outfile->Close();
    infile->Close();
}
