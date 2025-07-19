// Skimming program based on skim.cc by Giovanni Abbiendi, Emma Hess

#include <algorithm>
#include <array>
#include <cmath>      // for std::abs and std::fabs
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "FairMCEventHeader.h"
#include "MUonETrack.h"
#include "Station_hits.h"
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
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TVector3.h"
#include "releff.h"
#include <TStyle.h>


using namespace std;

// Compute ε = num/den and σ = sqrt(ε(1-ε)/den), then print label, ε±σ
void printEff(const char *label, double num, double den)
{
  if (den <= 0)
  {
    cout << label << ":  —  (denominator zero or negative)\n";
    return;
  }
  double eps = num / den;
  double err = sqrt(eps * (1.0 - eps) / den);
  cout
      << label
      << ":  " << num << "/" << den
      << "  =  " << eps
      << " ± " << err
      << "\n";
}

// added the whole tree
int main(int argc, char *argv[])
{

  int NMODULES = 18;
  string filelist = "runs/";
  string inputdir = "inputs/";
  string input;
  int iskim;
  string skimname;
  string outdir = "outputs/";

  if (argc != 4)
  {
    cout << "ERROR! Usage: " << argv[0] << " <run_nr>  <input (merged file list)>  <skim_code> " << endl;
    return 10;
  }
  else
  {
    filelist = filelist + argv[1] + "_files.txt";
    inputdir = inputdir + argv[1] + "/";
    input = argv[2];
    iskim = std::stoi(argv[3]);
    outdir += argv[1];

    if (iskim == 0)
    {
      cout << "\n SKIM: Analysis without saving skims \n"
           << endl;
    }
    else if (iskim == 1)
    {
      cout << "\n SKIM: single-mu interactions \n"
           << endl;
      skimname = "skim_interactions_singlemu";
    }
    else if (iskim == 60)
    {
      cout << "\n SKIM: pileup (2,3,4) mu interactions \n"
           << endl;
      skimname = "skim_interactions_pileup234mu";
    }
    else if (iskim == 61)
    {
      cout << "\n SKIM: single-mu OR pileup (2,3,4) mu interactions \n"
           << endl;
      skimname = "skim_interactions_any1234mu";
    }
    else if (iskim == 1001)
    {
      cout << "\n SKIM: single-mu interactions for MC \n"
           << endl;
      skimname = "skim_interactions_singlemu_MC";
    }
    else
    {
      cout << "\n NO OUTPUT SKIM selected! \n"
           << endl;
      return 20;
    }
  }

  //////////////////////////////////////////////////////////////////////
  // INPUT RUN full file list
  //
  cout << "\n *** INPUT full file list: " << filelist << endl;
  ifstream run_file_list(filelist);
  string i_file;
  vector<string> files;
  while (getline(run_file_list, i_file))
  {
    files.push_back(i_file);
    cout << i_file << endl;
  }
  cout << argv[1] << " filelist contains " << files.size() << " files " << endl;

  // INPUT RUN worked-on file list
  //
  TChain chain("cbmsim");
  //
  cout << "\n *** INPUT worked file list: " << inputdir + input << endl;
  ifstream input_file_list(inputdir + input);
  vector<string> input_files;
  while (getline(input_file_list, i_file))
  {
    input_files.push_back(i_file);
    string filepath = "/eos/user/p/pasenov/www/pasenov/MUonE_data/Condor_results/3stations_2025/MCminBias_3stations/";
    filepath = filepath + argv[1] + "/" + i_file;
    chain.Add(filepath.c_str());
    cout << i_file << endl;
    cout << filepath << endl;
  }

  cout << "chained " << input_files.size() << " files " << endl;
  auto n_entries = chain.GetEntries();
  cout << "Total number of input events  = " << n_entries << endl;

  // Declaration of leaf types
  vector<float> *Bend = 0;
  vector<unsigned short> *Bx = 0;
  vector<unsigned short> *Link = 0;
  vector<float> *LocalX = 0;
  vector<float> *LocalY = 0;
  vector<unsigned short> *StationID = 0;
  vector<unsigned int> *SuperID = 0;

  TClonesArray *MCTrack = nullptr;
  TClonesArray *TrackerStubs = nullptr;


  // Define a struct to hold skim flags
  struct OfflineSkimMode_t
  {
    Int_t StationIndex;
    Bool_t SingleMu;
    Bool_t Pileup234Mu;
  } OfflineSkimMode;

  chain.SetBranchAddress("Bend", &Bend);
  chain.SetBranchAddress("Bx", &Bx);
  chain.SetBranchAddress("Link", &Link);
  chain.SetBranchAddress("LocalX", &LocalX);
  chain.SetBranchAddress("LocalY", &LocalY);
  chain.SetBranchAddress("StationID", &StationID);

  chain.SetBranchAddress("MCTrack", &MCTrack);
  chain.SetBranchAddress("TrackerStubs", &TrackerStubs);

  TFile *o_file;
  TTree *o_tree = nullptr;
  Long64_t o_nevt = 0;

  if (iskim > 0 && iskim != 1001)
  {
    //
    // OUTPUT
    //
    // max size of skimmed output files: 2 GB
    TTree::SetMaxTreeSize(2000000000LL);
    o_file = TFile::Open((outdir + "/" + skimname + "_" + argv[1] + "_" + argv[2] + ".root").c_str(), "recreate");
    if (!o_file)
    {
      cout << "Error opening file" << std::endl;
      exit(-1);
    }

    // Clone the structure of the input tree, but with zero entries
    o_tree = chain.CloneTree(0);
    // o_tree = new TTree("cbmsim", "");

    o_tree->Branch("Bend", &Bend);
    o_tree->Branch("Bx", &Bx);
    o_tree->Branch("Link", &Link);
    o_tree->Branch("SuperID", &SuperID);
    // o_tree->Branch("UserBits", &o_UserBits);
    o_tree->Branch("LocalX", &LocalX);
    o_tree->Branch("LocalY", &LocalY);
    o_tree->Branch("StationID", &StationID);

    o_tree->Branch("MCTrack", &MCTrack);
    o_tree->Branch("TrackerStubs", &TrackerStubs);

    ///////// my new leafs
    /////////////////////////////////////////Maybe not/////////////////

    o_tree->Branch("OfflineSkimMode", &OfflineSkimMode, "StationIndex/I:SingleMu/O:Pileup234Mu/O");

    //////////////////////////////////////
    if (iskim == 1001)
    {
      TClonesArray *MCTrack = nullptr;
      TClonesArray *TrackerStubs = nullptr;
      // FairMCEventHeader *MCEventHeader = nullptr;

      chain.SetBranchAddress("MCTrack", &MCTrack);
      chain.SetBranchAddress("TrackerStubs", &TrackerStubs);
      // chain.SetBranchAddress("MCEventHeader.", &MCEventHeader);

      o_tree->Branch("MCTrack", &MCTrack);
      o_tree->Branch("TrackerStubs", &TrackerStubs);
      // o_tree->Branch("MCEventHeader.", &MCEventHeader);
    }
  }

  //////////////////////////////////////////////////////////////////////

  // histograms: find X range (file number - SuperId - time)
  //


  auto pos1 = files.at(0).find('-');
  auto pos2 = files.at(0).find('.', pos1); // trova il punto dopo il '-'

  unsigned int sid_min = std::stoul(files.at(0).substr(pos1 + 1, pos2 - pos1 - 1));
  unsigned int sid_max = sid_min;

  for (const auto &f : files)
  {


    pos1 = files.at(0).find('-');
    pos2 = files.at(0).find('.', pos1); // trova il punto dopo il '-'

    unsigned int sid_0 = std::stoul(f.substr(pos1 + 1, pos2 - pos1 - 1));

    sid_min = std::min(sid_min, sid_0);
    sid_max = std::max(sid_max, sid_0);
  }
  sid_max += 10000; // current decoding in chunks of 10000 SuperID's
  cout << "sid_min, sid_max = " << sid_min << ", " << sid_max << endl;
  int nbinsx = (sid_max - sid_min) / 10000;
  cout << "number of x bins = " << nbinsx << endl;

  // assume maximum 70 stubs per event in the two stations
  // (the real maximum in 2023 is 60, the DAQ cuts off hits beyond this limit)
  const Double_t nmax = 70;
  Int_t nb = nmax;

  string hname = outdir + "/histos_" + argv[1] + "_" + argv[2] + ".root";
  TFile *hfile = new TFile(hname.c_str(), "RECREATE");

  // stub counts per module and station
  TH1D *h_nstubsPerModule = new TH1D("h_nstubsPerModule", "nstubs on modules", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_S0 = new TH1D("h_nstubsPerModule_S0", "nstubs on modules - Station 0", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_S1 = new TH1D("h_nstubsPerModule_S1", "nstubs on modules - Station 1", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_S2 = new TH1D("h_nstubsPerModule_S2", "nstubs on modules - Station 2", NMODULES, 0, NMODULES);
  // for different selections
  TH1D *h_nstubsPerModule_presel = new TH1D("h_nstubsPerModule_presel", "nstubs on modules (preselected events)", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_trackable_0_onehit_1 = new TH1D("h_nstubsPerModule_trackable_0_onehit_1", "nstubs on modules (S0_trackable_S1_onehit events))", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_trackable_1_onehit_0 = new TH1D("h_nstubsPerModule_trackable_1_onehit_0", "nstubs on modules (S1_trackable_S0_onehit events))", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_trackable = new TH1D("h_nstubsPerModule_trackable", "nstubs on modules (trackable events)", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_fired12 = new TH1D("h_nstubsPerModule_fired12", "nstubs on modules (fired12 events)", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_single_clean = new TH1D("h_nstubsPerModule_single_clean", "nstubs on modules (clean single mu interaction events)", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_pileup234_skim = new TH1D("h_nstubsPerModule_pileup234_skim", "nstubs on modules (pileup234_skim events)", NMODULES, 0, NMODULES);
  TH1D *h_nstubsPerModule_fired12_plus_single_clean = new TH1D("h_nstubsPerModule_fired12_plus_single_clean", "nstubs on modules (fired12_plus_single_clean events)", NMODULES, 0, NMODULES);

  // for Luminosity: preselected single muon plots of Y % X position in station 0 (using last pair of XY modules)
  const Double_t nstrips = 1016;
  Int_t nbxy = nstrips / 8;
  TH2D *h2_LocalXY_presel_passing_golden = new TH2D("h2_LocalXY_presel_passing_golden", "LocalY vs LocalX (presel passing golden in S0)", nbxy, 0, nstrips, nbxy, 0, nstrips);

  // stub multiplicities for different selections
  TH1D *h_nTotStubs = new TH1D("h_nTotStubs", "Total N of Stubs; Number of Stubs; Events", nb, 0, nmax);
  TH2D *h2_nStubs = new TH2D("h2_nStubs", "N stubs on S1 vs S0", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0 = new TH1D("h_nStubs_0", "Tot N stubs First Station - S0", nb, 0, nmax);
  TH1D *h_nStubs_1 = new TH1D("h_nStubs_1", "Tot N stubs Second Station - S1", nb, 0, nmax);
  TH1D *h_nStubs_2 = new TH1D("h_nStubs_2", "Tot N stubs Second Station - S2", nb, 0, nmax);

  TH1D *h_nTotStubs_presel = new TH1D("h_nTotStubs_presel", "Total N of Stubs (preselected events); Number of Stubs; Events", nb, 0, nmax);
  TH2D *h2_nStubs_presel = new TH2D("h2_nStubs_presel", "N stubs on S1 vs S0 (preselected events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_presel = new TH1D("h_nStubs_0_presel", "Tot N stubs First Station - S0 (preselected events)", nb, 0, nmax);
  TH1D *h_nStubs_1_presel = new TH1D("h_nStubs_1_presel", "Tot N stubs Second Station - S1 (preselected events)", nb, 0, nmax);

  TH1D *h_nTotStubs_S0_trackable_S1_onehit = new TH1D("h_nTotStubs_S0_trackable_S1_onehit", "Total N of Stubs (S0_trackable_S1_onehit events", nb, 0, nmax);
  TH2D *h2_nStubs_S0_trackable_S1_onehit = new TH2D("h2_nStubs_S0_trackable_S1_onehit", "N stubs on S1 vs S0 (S0_trackable_S1_onehit events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_S0_trackable_S1_onehit = new TH1D("h_nStubs_0_S0_trackable_S1_onehit", "Tot N stubs First Station (S0_trackable_S1_onehit events)", nb, 0, nmax);
  TH1D *h_nStubs_1_S0_trackable_S1_onehit = new TH1D("h_nStubs_1_S0_trackable_S1_onehit", "Tot N stubs Second Station (S0_trackable_S1_onehit events)", nb, 0, nmax);

  TH1D *h_nTotStubs_S1_trackable_S0_onehit = new TH1D("h_nTotStubs_S1_trackable_S0_onehit", "Total N of Stubs (S1_trackable_S0_onehit events", nb, 0, nmax);
  TH2D *h2_nStubs_S1_trackable_S0_onehit = new TH2D("h2_nStubs_S1_trackable_S0_onehit", "N stubs on S1 vs S0 (S1_trackable_S0_onehit events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_S1_trackable_S0_onehit = new TH1D("h_nStubs_0_S1_trackable_S0_onehit", "Tot N stubs First Station (S1_trackable_S0_onehit events)", nb, 0, nmax);
  TH1D *h_nStubs_1_S1_trackable_S0_onehit = new TH1D("h_nStubs_1_S1_trackable_S0_onehit", "Tot N stubs Second Station (S1_trackable_S0_onehit events)", nb, 0, nmax);

  TH1D *h_nTotStubs_trackable = new TH1D("h_nTotStubs_trackable", "Total N of Stubs (trackable events); Number of Stubs; Events", nb, 0, nmax);
  TH2D *h2_nStubs_trackable = new TH2D("h2_nStubs_trackable", "N stubs on S1 vs S0 (trackable events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_trackable = new TH1D("h_nStubs_0_trackable", "Tot N stubs First Station - S0 (trackable events)", nb, 0, nmax);
  TH1D *h_nStubs_1_trackable = new TH1D("h_nStubs_1_trackable", "Tot N stubs Second Station - S1 (trackable events)", nb, 0, nmax);

  TH1D *h_nTotStubs_zero_0 = new TH1D("h_nTotStubs_zero_0", "Total N of Stubs - events with zero hits in S_0", nb, 0, nmax);
  TH1D *h_nTotStubs_zero_1 = new TH1D("h_nTotStubs_zero_1", "Total N of Stubs - events with zero hits in S_1", nb, 0, nmax);
  TH1D *h_nTotStubs_noise = new TH1D("h_nTotStubs_noise", "Total N of Stubs noise events", nb, 0, nmax);
  TH1D *h_nTotStubs_missing_0 = new TH1D("h_nTotStubs_missing_0", "Total N of Stubs muons missing S_0", nb, 0, nmax);
  TH1D *h_nTotStubs_missing_1 = new TH1D("h_nTotStubs_missing_1", "Total N of Stubs muons missing S_1", nb, 0, nmax);
  TH1D *h_nTotStubs_not_trackable = new TH1D("h_nTotStubs_not_trackable", "Total N of Stubs - not trackable events", nb, 0, nmax);
  TH1D *h_nTotStubs_presel_not_trackable = new TH1D("h_nTotStubs_presel_not_trackable", "Total N of Stubs (preselected events) - not trackable events", nb, 0, nmax);

  TH1D *h_nTotStubs_passing_1mu = new TH1D("h_nTotStubs_passing_1mu", "Total N of Stubs for passing mu", nb, 0, nmax);
  TH1D *h_nTotStubs_passing_2mu = new TH1D("h_nTotStubs_passing_2mu", "Total N of Stubs for passing 2mu", nb, 0, nmax);
  TH1D *h_nTotStubs_passing_3mu = new TH1D("h_nTotStubs_passing_3mu", "Total N of Stubs for passing 3mu", nb, 0, nmax);
  TH1D *h_nTotStubs_passing_4mu = new TH1D("h_nTotStubs_passing_4mu", "Total N of Stubs for passing 4mu", nb, 0, nmax);

  TH1D *h_nTotStubs_golden = new TH1D("h_nTotStubs_golden", "Total N of Stubs for Golden signal selection; Number of Stubs; Events", nb, 0, nmax);
  TH1D *h_nTotStubs_umsel = new TH1D("h_nTotStubs_umsel", "Total N of Stubs for Umberto's selection", nb, 0, nmax);
  TH1D *h_nTotStubs_single_cand = new TH1D("h_nTotStubs_single_cand", "Total N of Stubs for candidate single mu signal", nb, 0, nmax);
  TH1D *h_nTotStubs_single_clean = new TH1D("h_nTotStubs_single_clean", "Total N of Stubs for clean single mu signal", nb, 0, nmax);
  TH1D *h_nTotStubs_single_clean_2 = new TH1D("h_nTotStubs_single_clean_2", "Total N of Stubs for clean single mu signal second pair", nb, 0, nmax);
  TH2D *h2_nStubs_single_clean = new TH2D("h2_nStubs_single_clean", "N stubs on S1 vs S0 (single_clean events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_single_clean = new TH1D("h_nStubs_0_single_clean", "Tot N stubs First Station - S0 (single_clean events)", nb, 0, nmax);
  TH1D *h_nStubs_1_single_clean = new TH1D("h_nStubs_1_single_clean", "Tot N stubs Second Station - S1 (single_clean events)", nb, 0, nmax);

  TH1D *h_nTotStubs_pileup_any = new TH1D("h_nTotStubs_pileup_any", "Total N of Stubs for signal+(any)mu pileup", nb, 0, nmax);
  TH1D *h_nTotStubs_pileup_2mu = new TH1D("h_nTotStubs_pileup_2mu", "Total N of Stubs for 2mu pileup", nb, 0, nmax);
  TH1D *h_nTotStubs_pileup_3mu = new TH1D("h_nTotStubs_pileup_3mu", "Total N of Stubs for 3mu pileup", nb, 0, nmax);
  TH1D *h_nTotStubs_pileup_4mu = new TH1D("h_nTotStubs_pileup_4mu", "Total N of Stubs for 4mu pileup", nb, 0, nmax);
  TH1D *h_nTotStubs_pileup_many = new TH1D("h_nTotStubs_pileup_many", "Total N of Stubs for >=4mu pileup", nb, 0, nmax);
  TH1D *h_nTotStubs_pileup_skim = new TH1D("h_nTotStubs_pileup_skim", "Total N of Stubs (pileup_skim events); Number of Stubs; Events", nb, 0, nmax);
  TH2D *h2_nStubs_pileup_skim = new TH2D("h2_nStubs_pileup_skim", "N stubs on S1 vs S0 (pileup_skim events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_pileup_skim = new TH1D("h_nStubs_0_pileup_skim", "Tot N stubs First Station - S0 (pileup_skim events)", nb, 0, nmax);
  TH1D *h_nStubs_1_pileup_skim = new TH1D("h_nStubs_1_pileup_skim", "Tot N stubs Second Station - S1 (pileup_skim events)", nb, 0, nmax);

  TH1D *h_nTotStubs_pileup234_skim = new TH1D("h_nTotStubs_pileup234_skim", "Total N of Stubs (pileup234_skim events); Number of Stubs; Events", nb, 0, nmax);
  TH2D *h2_nStubs_pileup234_skim = new TH2D("h2_nStubs_pileup234_skim", "N stubs on S1 vs S0 (pileup234_skim events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_pileup234_skim = new TH1D("h_nStubs_0_pileup234_skim", "Tot N stubs First Station - S0 (pileup234_skim events)", nb, 0, nmax);
  TH1D *h_nStubs_1_pileup234_skim = new TH1D("h_nStubs_1_pileup234_skim", "Tot N stubs Second Station - S1 (pileup234_skim events)", nb, 0, nmax);

  TH1D *h_nTotStubs_fired12 = new TH1D("h_nTotStubs_fired12", "Total N of Stubs (fired12 events); Number of Stubs; Events", nb, 0, nmax);
  TH2D *h2_nStubs_fired12 = new TH2D("h2_nStubs_fired12", "N stubs on S1 vs S0 (fired12 events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_fired12 = new TH1D("h_nStubs_0_fired12", "Tot N stubs First Station - S0 (fired12 events)", nb, 0, nmax);
  TH1D *h_nStubs_1_fired12 = new TH1D("h_nStubs_1_fired12", "Tot N stubs Second Station - S1 (fired12 events)", nb, 0, nmax);

  TH1D *h_nTotStubs_fired12_plus_single_clean = new TH1D("h_nTotStubs_fired12_plus_single_clean", "Total N of Stubs (fired12_plus_single_clean events); Number of Stubs; Events", nb, 0, nmax);
  TH2D *h2_nStubs_fired12_plus_single_clean = new TH2D("h2_nStubs_fired12_plus_single_clean", "N stubs on S1 vs S0 (fired12_plus_single_clean events)", nb, 0, nmax, nb, 0, nmax);
  TH1D *h_nStubs_0_fired12_plus_single_clean = new TH1D("h_nStubs_0_fired12_plus_single_clean", "Tot N stubs First Station - S0 (fired12_plus_single_clean events)", nb, 0, nmax);
  TH1D *h_nStubs_1_fired12_plus_single_clean = new TH1D("h_nStubs_1_fired12_plus_single_clean", "Tot N stubs Second Station - S1 (fired12_plus_single_clean events)", nb, 0, nmax);

  // FOR STATION EFFICIENCIES (efficiency of Trackable Track pattern)
  // time profile histos for all recorded events
  TH1D *h_goldenS0_vs_time = new TH1D("h_goldenS0_vs_time", "# golden S_0 events vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_goldenS0_trackableS1_vs_time = new TH1D("h_goldenS0_trackableS1_vs_time", "# golden S_0 events trackable in S_1 vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_effS1_vs_time = new TH1D("h_effS1_vs_time", "S_1 Efficiency vs time", nbinsx, 0., double(nbinsx));
  //
  TH1D *h_goldenS1_vs_time = new TH1D("h_goldenS1_vs_time", "# golden S_1 events vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_goldenS1_trackableS0_vs_time = new TH1D("h_goldenS1_trackableS0_vs_time", "# golden S_1 events trackable in S_0 vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_effS0_vs_time = new TH1D("h_effS0_vs_time", "S_0 Efficiency vs time", nbinsx, 0., double(nbinsx));
  //
  // time profile histos for only preselected events (minimum one hit requested in the station under test)
  TH1D *h_presel_goldenS0_vs_time = new TH1D("h_presel_goldenS0_vs_time", "# golden S_0 events vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_effS1_vs_time = new TH1D("h_presel_effS1_vs_time", "S_1 Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  //
  TH1D *h_presel_goldenS1_vs_time = new TH1D("h_presel_goldenS1_vs_time", "# golden S_1 events vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_effS0_vs_time = new TH1D("h_presel_effS0_vs_time", "S_0 Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  //
  // (efficiency of Golden Track pattern)
  TH1D *h_golden_effS1_vs_time = new TH1D("h_golden_effS1_vs_time", "S_1 Golden Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_golden_effS0_vs_time = new TH1D("h_golden_effS0_vs_time", "S_0 Golden Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  //
  // ******************* 2 MU efficiencies *******************
  // ******************* efficiencies in S1
  TH1D *h_presel_passing_2mu_goldenS0_vs_time = new TH1D("h_presel_passing_2mu_goldenS0_vs_time", "# golden 2mu S_0 events vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_passing_2mu_goldenS0_trackableS1_vs_time = new TH1D("h_presel_passing_2mu_goldenS0_trackableS1_vs_time", "# golden 2mu S_0 events trackable in S_1 vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_2mu_effS1_vs_time = new TH1D("h_presel_2mu_effS1_vs_time", "S_1 2mu Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_golden_2mu_effS1_vs_time = new TH1D("h_golden_2mu_effS1_vs_time", "S_1 2mu Golden Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  // ******************* efficiencies in S0
  TH1D *h_presel_passing_2mu_goldenS1_vs_time = new TH1D("h_presel_passing_2mu_goldenS1_vs_time", "# golden 2mu S_1 events vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_passing_2mu_goldenS1_trackableS0_vs_time = new TH1D("h_presel_passing_2mu_goldenS1_trackableS0_vs_time", "# golden 2mu S_1 events trackable in S_0 vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_2mu_effS0_vs_time = new TH1D("h_presel_2mu_effS0_vs_time", "S_0 2mu Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_golden_2mu_effS0_vs_time = new TH1D("h_golden_2mu_effS0_vs_time", "S_0 2mu Golden Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  // ******************* 3 MU efficiencies *******************
  // ******************* efficiencies in S1
  TH1D *h_presel_passing_3mu_goldenS0_vs_time = new TH1D("h_presel_passing_3mu_goldenS0_vs_time", "# golden 3mu S_0 events vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_passing_3mu_goldenS0_trackableS1_vs_time = new TH1D("h_presel_passing_3mu_goldenS0_trackableS1_vs_time", "# golden 3mu S_0 events trackable in S_1 vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_3mu_effS1_vs_time = new TH1D("h_presel_3mu_effS1_vs_time", "S_1 3mu Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_golden_3mu_effS1_vs_time = new TH1D("h_golden_3mu_effS1_vs_time", "S_1 3mu Golden Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  // ******************* efficiencies in S0
  TH1D *h_presel_passing_3mu_goldenS1_vs_time = new TH1D("h_presel_passing_3mu_goldenS1_vs_time", "# golden 3mu S_1 events vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_passing_3mu_goldenS1_trackableS0_vs_time = new TH1D("h_presel_passing_3mu_goldenS1_trackableS0_vs_time", "# golden 3mu S_1 events trackable in S_0 vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_3mu_effS0_vs_time = new TH1D("h_presel_3mu_effS0_vs_time", "S_0 3mu Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_golden_3mu_effS0_vs_time = new TH1D("h_golden_3mu_effS0_vs_time", "S_0 3mu Golden Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  // ******************* 4 MU efficiencies *******************
  // ******************* efficiencies in S1
  TH1D *h_presel_passing_4mu_goldenS0_vs_time = new TH1D("h_presel_passing_4mu_goldenS0_vs_time", "# golden 4mu S_0 events vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_passing_4mu_goldenS0_trackableS1_vs_time = new TH1D("h_presel_passing_4mu_goldenS0_trackableS1_vs_time", "# golden 4mu S_0 events trackable in S_1 vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_4mu_effS1_vs_time = new TH1D("h_presel_4mu_effS1_vs_time", "S_1 4mu Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_golden_4mu_effS1_vs_time = new TH1D("h_golden_4mu_effS1_vs_time", "S_1 4mu Golden Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  // ******************* efficiencies in S0
  TH1D *h_presel_passing_4mu_goldenS1_vs_time = new TH1D("h_presel_passing_4mu_goldenS1_vs_time", "# golden 4mu S_1 events vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_passing_4mu_goldenS1_trackableS0_vs_time = new TH1D("h_presel_passing_4mu_goldenS1_trackableS0_vs_time", "# golden 4mu S_1 events trackable in S_0 vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_presel_4mu_effS0_vs_time = new TH1D("h_presel_4mu_effS0_vs_time", "S_0 4mu Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));
  TH1D *h_golden_4mu_effS0_vs_time = new TH1D("h_golden_4mu_effS0_vs_time", "S_0 4mu Golden Efficiency vs time (after preselection)", nbinsx, 0., double(nbinsx));

  // event RATE plots for categories
  //
  TH1D *h_passing_1mu_golden_vs_time = new TH1D("h_passing_1mu_golden_vs_time", "# passing 1mu golden events vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_passing_2mu_golden_vs_time = new TH1D("h_passing_2mu_golden_vs_time", "# passing 2mu golden events vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_passing_3mu_golden_vs_time = new TH1D("h_passing_3mu_golden_vs_time", "# passing 3mu golden events vs time", nbinsx, 0., double(nbinsx));
  TH1D *h_passing_4mu_golden_vs_time = new TH1D("h_passing_4mu_golden_vs_time", "# passing 4mu golden events vs time", nbinsx, 0., double(nbinsx));
  
  Long64_t N_input_events = 0;
  Long64_t N_input_stubs_0 = 0;
  Long64_t N_input_stubs_1 = 0;
  Long64_t N_input_stubs_2 = 0;
  Long64_t N_input_stubs = 0; // total S0+S1
  Long64_t N_stubs_skim = 0;

  Int_t N_single_muon_offline_online_0 = 0;
  Int_t N_single_muon_offline_online_1 = 0;

  Int_t N_pileup_muon_offline_online_0 = 0;
  Int_t N_pileup_muon_offline_online_1 = 0;

  Int_t N_single_passing_muon_offline_online_0 = 0;
  Int_t N_pileup_passing_muon_offline_online_0 = 0;
  Int_t N_single_passing_muon_offline_online_1 = 0;
  Int_t N_pileup_passing_muon_offline_online_1 = 0;

  Long64_t N_presel_events = 0;
  Long64_t N_presel_stubs_0 = 0;
  Long64_t N_presel_stubs_1 = 0;
  Long64_t N_presel_stubs = 0; // total S0+S1

  //  Long64_t N_trackable_0_onehit_1 = 0;
  Long64_t N_trackable_0_onehit_1_stubs_0 = 0;
  Long64_t N_trackable_0_onehit_1_stubs_1 = 0;
  Long64_t N_trackable_0_onehit_1_stubs = 0; // total S0+S1

  //  Long64_t N_trackable_1_onehit_0 = 0;
  Long64_t N_trackable_1_onehit_0_stubs_0 = 0;
  Long64_t N_trackable_1_onehit_0_stubs_1 = 0;
  Long64_t N_trackable_1_onehit_0_stubs = 0; // total S0+S1

  //  Long64_t N_single_clean_events = 0;
  Long64_t N_single_clean_stubs_0 = 0;
  Long64_t N_single_clean_stubs_1 = 0;
  Long64_t N_single_clean_stubs = 0; // total S0+S1

  //  Long64_t N_pileup_skim_events = 0;
  Long64_t N_pileup_skim_stubs_0 = 0;
  Long64_t N_pileup_skim_stubs_1 = 0;
  Long64_t N_pileup_skim_stubs = 0; // total S0+S1

  //  Long64_t N_pileup234_skim_events = 0;
  Long64_t N_pileup234_skim_stubs_0 = 0;
  Long64_t N_pileup234_skim_stubs_1 = 0;
  Long64_t N_pileup234_skim_stubs = 0; // total S0+S1

  Long64_t N_zero_0 = 0;
  Long64_t N_noise_0 = 0;
  //
  Long64_t N_onehit_0 = 0;
  Long64_t N_twohit_0 = 0;
  Long64_t N_threehit_0 = 0;
  Long64_t N_fourhit_0 = 0;
  Long64_t N_trackable_0 = 0;
  Long64_t N_trackable_0_onehit_1 = 0;
  Long64_t N_trackable_0_twohit_1 = 0;
  Long64_t N_trackable_0_threehit_1 = 0;
  Long64_t N_trackable_0_fourhit_1 = 0;
  Long64_t N_trackable_1_onehit_0 = 0;
  Long64_t N_trackable_events = 0;
  //
  Long64_t N_passing_1mu_golden_0 = 0;
  Long64_t N_passing_2mu_golden_0 = 0;
  Long64_t N_passing_3mu_golden_0 = 0;
  Long64_t N_passing_4mu_golden_0 = 0;
  Long64_t N_single_cand_0 = 0;
  Long64_t N_single_cand_1 = 0;
  Long64_t N_single_clean_0 = 0;
  Long64_t N_single_clean_1 = 0;
  Long64_t N_pileup_any_0 = 0;
  Long64_t N_pileup_2mu_0 = 0;
  Long64_t N_pileup_3mu_0 = 0;
  Long64_t N_pileup_4mu_0 = 0;
  Long64_t N_pileup_many_0 = 0;

  Long64_t N_zero_1 = 0;
  Long64_t N_noise_1 = 0;
  Long64_t N_passing_1mu_golden_1 = 0;
  Long64_t N_passing_2mu_golden_1 = 0;
  Long64_t N_passing_3mu_golden_1 = 0;
  Long64_t N_passing_4mu_golden_1 = 0;
  Long64_t N_passing_cand_1 = 0;
  Long64_t N_passing_clean_1 = 0;
  Long64_t N_2tracks_1 = 0;
  Long64_t N_2tracks_2 = 0;
  Long64_t N_golden_1 = 0;

  Long64_t N_signal_plus_pileup_1 = 0;
  Long64_t N_signal_plus_tripileup_1 = 0;
  Long64_t N_signal_plus_fourpileup_1 = 0;
  Long64_t N_pileup_2mu_1 = 0;
  Long64_t N_pileup_3mu_1 = 0;
  Long64_t N_pileup_4mu_1 = 0;

  Long64_t N_zero_events = 0;
  Long64_t N_noise_events = 0;
  Long64_t N_missing_0_events = 0;
  Long64_t N_missing_1_events = 0;
  Long64_t N_presel_not_trackable_events = 0;
  Long64_t N_not_trackable_events = 0;
  Long64_t N_passing_golden_S0_zero_S1_events = 0;
  Long64_t N_passing_golden_S1_zero_S0_events = 0;

  Long64_t N_passing_1mu_events = 0;
  Long64_t N_passing_2mu_events = 0;
  Long64_t N_passing_3mu_events = 0;
  Long64_t N_passing_4mu_events = 0;

  Long64_t N_passing_1mu_golden_events = 0;
  Long64_t N_passing_2mu_golden_events = 0;
  Long64_t N_passing_3mu_golden_events = 0;
  Long64_t N_passing_4mu_golden_events = 0;

  Long64_t N_passing_golden_S0_trackable_S1_events = 0;
  Long64_t N_passing_golden_S1_trackable_S0_events = 0;
  Long64_t N_presel_passing_golden_S0_events = 0;
  Long64_t N_presel_passing_golden_S1_events = 0;
  //
  Long64_t N_presel_passing_2mu_golden_S0_events = 0;
  Long64_t N_presel_passing_2mu_golden_S1_events = 0;
  Long64_t N_presel_passing_2mu_golden_S0_trackable_S1_events = 0;
  Long64_t N_presel_passing_2mu_golden_S1_trackable_S0_events = 0;
  //
  Long64_t N_presel_passing_3mu_golden_S0_events = 0;
  Long64_t N_presel_passing_3mu_golden_S1_events = 0;
  Long64_t N_presel_passing_3mu_golden_S0_trackable_S1_events = 0;
  Long64_t N_presel_passing_3mu_golden_S1_trackable_S0_events = 0;
  //
  Long64_t N_presel_passing_4mu_golden_S0_events = 0;
  Long64_t N_presel_passing_4mu_golden_S1_events = 0;
  Long64_t N_presel_passing_4mu_golden_S0_trackable_S1_events = 0;
  Long64_t N_presel_passing_4mu_golden_S1_trackable_S0_events = 0;

  Long64_t N_golden_events = 0;
  Long64_t N_umsel_events = 0;

  Long64_t N_single_cand_events = 0;
  Long64_t N_single_clean_events = 0;
  Long64_t N_single_clean_events_1 = 0;
  Long64_t N_pileup_any_events = 0;
  Long64_t N_pileup_2mu_events = 0;
  Long64_t N_pileup_3mu_events = 0;
  Long64_t N_pileup_4mu_events = 0;
  Long64_t N_pileup_many_events = 0;
  Long64_t N_pileup_skim_events = 0;
  Long64_t N_pileup234_skim_events = 0;
  Long64_t N_loose_events = 0;
  Long64_t N_loose234_events = 0;

  Long64_t N_fired6_0 = 0;
  Long64_t N_fired6_1 = 0;
  Long64_t N_fired12_events = 0;
  Long64_t N_fired12_plus_single_clean_events = 0;


  Int_t N_signal_plus_pileup_2 = 0;
  Int_t N_pileup_2mu_events_1 = 0;
  Int_t N_signal_plus_tripileup_2 = 0;
  // Int_t N_pileup_3mu_2 = 0;
  Int_t N_pileup_3mu_events_1 = 0;
  Int_t N_signal_plus_fourpileup_2 = 0;
  // Int_t N_pileup_4mu_2 = 0;
  Int_t N_pileup_4mu_events_1 = 0;
  Int_t N_pileup234_skim_events_1 = 0;

  Long64_t N_passing_1mu_golden_events_1 = 0;
  Int_t N_passing_1mu_golden_2 = 0;
  Int_t N_passing_2mu_golden_events_1 = 0;
  Int_t N_passing_3mu_golden_events_1 = 0;
  Int_t N_passing_4mu_golden_events_1 = 0;
  Int_t N_passing_2mu_golden_2 = 0;
  Int_t N_passing_3mu_golden_2 = 0;
  Int_t N_passing_4mu_golden_2 = 0;
  Int_t N_empty = 0;
  Int_t N_differentSize = 0;

  for (Long64_t i = 0; i < n_entries; i++) // n_entries
  {

    Long64_t ientry = chain.LoadTree(i);
    if (ientry < 0)
    {
      cerr << "***ERROR in reading input chain, event = " << i << ", ientry = " << ientry << endl;
      break;
    }
    if (i % 1000 == 0)
      cout << " processing event : " << i << "\r" << flush << endl;

    chain.GetEntry(i);
    N_input_events++;

    OfflineSkimMode.StationIndex = -1;
    OfflineSkimMode.SingleMu = false;
    OfflineSkimMode.Pileup234Mu = false;
    // time from beginning of run expressed in units of 10^4 orbits ~ 3564*25ns*10000 = 0.891 s
    // in this way each file should correspond to a bin in histograms

    int nTotStubs = Link->size(); // vectors size

    N_input_stubs += nTotStubs;

    h_nTotStubs->Fill(nTotStubs);


    std::array<int, 18> nstubs = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int nStubs_0 = 0;
    int nStubs_1 = 0;
    int nStubs_2 = 0;

    // size_t expected_size = tracker_Link->size();
    // Check if sizes match exactly
    size_t tracker_size = Bend->size();
    if (Bend->empty() || Link->empty() || LocalX->empty() || LocalY->empty() || StationID->empty())
    {
      cout << "Skipping event " << i << ": One or more tracker vectors are empty." << endl;
      if (Link->size() != tracker_size || LocalX->size() != tracker_size ||
          LocalY->size() != tracker_size || StationID->size() != tracker_size //||tracker_BxOffset->size() != tracker_size
      )
      {
        std::cout << "This event also has different sizes of ntuples " << std::endl;
      }

      //cout << "tracker_Bend size: " << tracker_Bend->size() << endl;
      N_empty++;
      continue;
    }

    ///////////////////////////////////////////////////////
    // loop on all stubs in the two stations

    for (int j = 0; j < nTotStubs; ++j)
    {
      // imod=(0:5) in station 0; imod=(6:11) in station 1
      int imod = Link->at(j);
      // int imod = Link->at(j) + StationID->at(j) * 6;

      if (imod < 0 || imod >= int(nstubs.size())) {
        std::cerr << "Warning: imod=" << imod << " out of bounds for nstubs (size=" << nstubs.size() << ") at event " << i << ", stub " << j << std::endl;
        continue;
      }

      nstubs.at(imod)++;

      h_nstubsPerModule->Fill(imod);

      if (imod < 6)
      {
        nStubs_0++;
        h_nstubsPerModule_S0->Fill(imod);
      }
      else if (imod < 12)
      {
        nStubs_1++;
        h_nstubsPerModule_S1->Fill(imod);
      }
      else
      {
        nStubs_2++;
        h_nstubsPerModule_S2->Fill(imod);
      }
    }

    // end loop on all stubs in the two stations
    //////////////////////// ///////////////////////////////

    N_input_stubs_0 += nStubs_0;
    N_input_stubs_1 += nStubs_1;
    N_input_stubs_2 += nStubs_2;

    h_nStubs_0->Fill(nStubs_0);
    h_nStubs_1->Fill(nStubs_1);
    h_nStubs_2->Fill(nStubs_2);
    h2_nStubs->Fill(nStubs_0, nStubs_1);

    //////////////////////////////////////////////////////////////
    // DEFINE the structures specifying the stations' hit patterns
    //////////////////////////////////////////////////////////////
    Station_hits S_0 = fill_station_hits(nstubs, 0);
    // cout << "Event " << i << ": S_0.multifired_xy_modules = " << S_0.multifired_xy_modules << endl;
    // cout << "Event " << i << ": S_0.multifired_uv_modules = " << S_0.multifired_uv_modules << endl;


    Station_hits S_1 = fill_station_hits(nstubs, 6);
    // cout << "Event " << i << ": S_1.multifired_xy_modules = " << S_1.multifired_xy_modules << endl;
    // cout << "Event " << i << ": S_1.multifired_uv_modules = " << S_1.multifired_uv_modules << endl;

    Station_hits S_2 = fill_station_hits(nstubs, 12);
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    bool ok_skim = false;
    bool ok_skim_MC = false;

    // conditions on station 0
    bool is_zero_0 = nStubs_0 == 0;
    if (is_zero_0)
    {
      N_zero_0++;
      h_nTotStubs_zero_0->Fill(nTotStubs);
    }

    bool is_noise_0 = (S_0.fired_xy_modules + S_0.fired_uv_modules) <= 1;
    if (is_noise_0)
      N_noise_0++;

    int fired_modules_0 = S_0.fired_xy_modules + S_0.fired_uv_modules;
    bool is_fired6_0 = fired_modules_0 == 6;
    if (is_fired6_0)
      N_fired6_0++;

    if (fired_modules_0 >= 1)
    {
      N_onehit_0++;
    }
    if (fired_modules_0 >= 2)
    {
      N_twohit_0++;
    }
    if (fired_modules_0 >= 3)
    {
      N_threehit_0++;
    }
    if (fired_modules_0 >= 4)
    {
      N_fourhit_0++;
    }

    bool is_trackable_0 = S_0.fired_xy_modules == 4 && S_0.fired_uv_modules > 0;
    if (is_trackable_0)
    {
      N_trackable_0++;
    }

    bool is_passing_golden_0 = nStubs_0 == 6 && S_0.fired_xy_modules == 4 && S_0.fired_uv_modules == 2;
    if (is_passing_golden_0)
    {

      N_passing_1mu_golden_0++;

      if (iskim == 100)
        ok_skim = true;
    }

    bool is_passing_2mu_golden_0 = nStubs_0 == 12 && S_0.multifired_xy_modules == 4 && S_0.multifired_uv_modules == 2;
    if (is_passing_2mu_golden_0)
    {
      N_passing_2mu_golden_0++;
    }

    bool is_passing_3mu_golden_0 = nStubs_0 == 18 && S_0.threefired_xy_modules == 4 && S_0.threefired_uv_modules == 2;
    if (is_passing_3mu_golden_0)
    {
      N_passing_3mu_golden_0++;
    }
    bool is_passing_4mu_golden_0 = nStubs_0 == 24 && S_0.fourfired_xy_modules == 4 && S_0.fourfired_uv_modules == 2;
    if (is_passing_4mu_golden_0)
    {
      N_passing_4mu_golden_0++;
    }

    bool is_passing_2mu_trackable_0 = S_0.multifired_xy_modules == 4 && S_0.nhits_uv >= 2;
    bool is_passing_3mu_trackable_0 = S_0.threefired_xy_modules == 4 && S_0.nhits_uv >= 3;
    bool is_passing_4mu_trackable_0 = S_0.fourfired_xy_modules == 4 && S_0.nhits_uv >= 4;

    bool is_single_cand_0 = is_trackable_0 && (S_0.multifired_xy_modules + S_0.multifired_uv_modules) <= 1;
    if (is_single_cand_0)
      N_single_cand_0++;
    bool is_single_clean_0 = is_single_cand_0 && nStubs_0 < 8;
    if (is_single_clean_0)
      N_single_clean_0++;

    bool is_fired6_plus_single_clean_0 =
        is_fired6_0 && S_0.multifired_xy_modules == 4 && S_0.multifired_uv_modules >= 1 && nStubs_0 <= 13;

    bool is_pileup_any_0 = S_0.multifired_xy_modules == 4 && S_0.nhits_uv >= 2;
    if (is_pileup_any_0)
      N_pileup_any_0++;

    bool is_pileup_2mu_0 = S_0.multifired_xy_modules == 4 && S_0.nhits_uv >= 2 &&
                           (S_0.threefired_xy_modules + S_0.threefired_uv_modules) <= 1 &&
                           (S_0.nhits_xy + S_0.nhits_uv) <= 13;
    if (is_pileup_2mu_0)
      N_pileup_2mu_0++;

    bool is_pileup_3mu_0 = S_0.threefired_xy_modules == 4 && S_0.nhits_uv >= 3 &&
                           (S_0.fourfired_xy_modules + S_0.fourfired_uv_modules) <= 1 &&
                           (S_0.nhits_xy + S_0.nhits_uv) <= 19;
    if (is_pileup_3mu_0)
      N_pileup_3mu_0++;

    bool is_pileup_4mu_0 = S_0.fourfired_xy_modules == 4 && S_0.nhits_uv >= 4 &&
                           (S_0.fivefired_xy_modules + S_0.fivefired_uv_modules) <= 1 &&
                           (S_0.nhits_xy + S_0.nhits_uv) <= 25;
    if (is_pileup_4mu_0)
      N_pileup_4mu_0++;

    // (4 or more muons)
    bool is_pileup_many_0 = S_0.fourfired_xy_modules == 4 && S_0.nhits_uv >= 4; // (nhits >=20)
    if (is_pileup_many_0)
      N_pileup_many_0++;

    // conditions on station 1
    //
    bool is_zero_1 = nStubs_1 == 0;
    if (is_zero_1)
    {
      N_zero_1++;
      h_nTotStubs_zero_1->Fill(nTotStubs);
    }

    bool is_noise_1 = (S_1.fired_xy_modules + S_1.fired_uv_modules) <= 1;
    if (is_noise_1)
      N_noise_1++;

    int fired_modules_1 = S_1.fired_xy_modules + S_1.fired_uv_modules;
    bool is_fired6_1 = fired_modules_1 == 6;
    if (is_fired6_1)
      N_fired6_1++;

    bool is_passing_golden_1 = nStubs_1 == 6 && S_1.fired_xy_modules == 4 && S_1.fired_uv_modules == 2;
    if (is_passing_golden_1)
    {
      N_passing_1mu_golden_1++;
    }
    bool is_passing_2mu_golden_1 = nStubs_1 == 12 && S_1.multifired_xy_modules == 4 && S_1.multifired_uv_modules == 2;
    if (is_passing_2mu_golden_1)
    {
      N_passing_2mu_golden_1++;
    }
    bool is_passing_3mu_golden_1 = nStubs_1 == 18 && S_1.threefired_xy_modules == 4 && S_1.threefired_uv_modules == 2;
    if (is_passing_3mu_golden_1)
    {
      N_passing_3mu_golden_1++;
    }
    bool is_passing_4mu_golden_1 = nStubs_1 == 24 && S_1.fourfired_xy_modules == 4 && S_1.fourfired_uv_modules == 2;
    if (is_passing_4mu_golden_1)
    {
      N_passing_4mu_golden_1++;
    }

    bool is_passing_2mu_golden_2 = nStubs_2 == 12 && S_2.multifired_xy_modules == 4 && S_2.multifired_uv_modules == 2;
    if (is_passing_2mu_golden_2)
    {
      N_passing_2mu_golden_2++;
    }
    bool is_passing_3mu_golden_2 = nStubs_2 == 18 && S_2.threefired_xy_modules == 4 && S_2.threefired_uv_modules == 2;
    if (is_passing_3mu_golden_2)
    {
      N_passing_3mu_golden_2++;
    }
    bool is_passing_4mu_golden_2 = nStubs_2 == 24 && S_2.fourfired_xy_modules == 4 && S_2.fourfired_uv_modules == 2;
    if (is_passing_4mu_golden_2)
    {
      N_passing_4mu_golden_2++;
    }

    bool is_passing_2mu_trackable_1 = S_1.multifired_xy_modules == 4 && S_1.nhits_uv >= 2;
    bool is_passing_3mu_trackable_1 = S_1.threefired_xy_modules == 4 && S_1.nhits_uv >= 3;
    bool is_passing_4mu_trackable_1 = S_1.fourfired_xy_modules == 4 && S_1.nhits_uv >= 4;

    bool is_trackable_1 = S_1.fired_xy_modules == 4 && S_1.fired_uv_modules > 0;
    bool is_passing_cand_1 = is_trackable_1 && (S_1.multifired_xy_modules + S_1.multifired_uv_modules) <= 1;
    if (is_passing_cand_1)
      N_passing_cand_1++;
    bool is_passing_clean_1 = is_passing_cand_1 && nStubs_1 < 8;
    if (is_passing_clean_1)
      N_passing_clean_1++;

    // second pair of stations, stat1
    bool is_single_cand_1 = is_trackable_1 && (S_1.multifired_xy_modules + S_1.multifired_uv_modules) <= 1;
    if (is_single_cand_1)
      N_single_cand_1++;

    bool is_single_clean_1 = is_single_cand_1 && nStubs_1 < 8;
    if (is_single_clean_1)
      N_single_clean_1++;

    bool is_2nd_pattern_1 = (S_1.multifired_xy_modules > 1 && S_1.nhits_uv > 1) ||
                            (S_1.multifired_xy_modules > 2 && S_1.nhits_uv > 0);
    bool is_2tracks_1 = is_trackable_1 && is_2nd_pattern_1;
    if (is_2tracks_1)
      N_2tracks_1++;

    bool is_fired6_plus_trackable_1 = is_fired6_1 && S_1.multifired_xy_modules == 4 && S_1.nhits_uv >= 3;

    bool is_2nd_pattern_overfired6_1 = (S_1.threefired_xy_modules >= 2 && S_1.nhits_uv >= 4) ||
                                       (S_1.threefired_xy_modules >= 3 && S_1.nhits_uv >= 3);

    bool is_fired6_2tracks_1 = is_fired6_plus_trackable_1 && is_2nd_pattern_overfired6_1;

    bool is_golden_1 = nStubs_1 == 12 && S_1.multifired_xy_modules == 4 && S_1.multifired_uv_modules == 2;
    if (is_golden_1)
      N_golden_1++;

    bool is_pileup_2mu_1 = S_1.multifired_xy_modules == 4 && S_1.nhits_uv >= 2 &&
                           (S_1.threefired_xy_modules + S_1.threefired_uv_modules) <= 1 &&
                           (S_1.nhits_xy + S_1.nhits_uv) <= 13;
    if (is_pileup_2mu_1)
      N_pileup_2mu_1++;

    bool is_pileup_3mu_1 = S_1.threefired_xy_modules == 4 && S_1.nhits_uv >= 3 &&
                           (S_1.fourfired_xy_modules + S_1.fourfired_uv_modules) <= 1 &&
                           (S_1.nhits_xy + S_1.nhits_uv) <= 19;
    if (is_pileup_3mu_1)
      N_pileup_3mu_1++;

    bool is_pileup_4mu_1 = S_1.fourfired_xy_modules == 4 && S_1.nhits_uv >= 4 &&
                           (S_1.fivefired_xy_modules + S_1.fivefired_uv_modules) <= 1 &&
                           (S_1.nhits_xy + S_1.nhits_uv) <= 25;
    if (is_pileup_4mu_1)
      N_pileup_4mu_1++;

    bool is_signal_plus_pileup_1 =
        S_1.multifired_xy_modules == 4 &&
        ((S_1.threefired_xy_modules >= 2 && S_1.nhits_uv >= 3) ||
         (S_1.threefired_xy_modules >= 3 && S_1.nhits_uv >= 2));
    if (is_signal_plus_pileup_1)
      N_signal_plus_pileup_1++;

    bool is_signal_plus_tripileup_1 =
        S_1.threefired_xy_modules == 4 &&
        ((S_1.fourfired_xy_modules >= 2 && S_1.nhits_uv >= 4) ||
         (S_1.fourfired_xy_modules >= 3 && S_1.nhits_uv >= 3));
    if (is_signal_plus_tripileup_1)
      N_signal_plus_tripileup_1++;

    bool is_signal_plus_fourpileup_1 =
        S_1.fourfired_xy_modules == 4 &&
        ((S_1.fivefired_xy_modules >= 2 && S_1.nhits_uv >= 5) ||
         (S_1.fivefired_xy_modules >= 3 && S_1.nhits_uv >= 4));
    if (is_signal_plus_fourpileup_1)
      N_signal_plus_fourpileup_1++;

    //////Third station

    bool is_trackable_2 = S_2.fired_xy_modules == 4 && S_2.fired_uv_modules > 0;

    bool is_2nd_pattern_2 = (S_2.multifired_xy_modules > 1 && S_2.nhits_uv > 1) ||
                            (S_2.multifired_xy_modules > 2 && S_2.nhits_uv > 0);

    bool is_2tracks_2 = is_trackable_2 && is_2nd_pattern_2;
    if (is_2tracks_2)
      N_2tracks_2++;

    // condition on S0+S1
    //
    bool is_zero_event = is_zero_0 && is_zero_1;
    if (is_zero_event)
      N_zero_events++;

    // at least one hit on both stations
    bool is_presel_event = !is_zero_0 && !is_zero_1;
    if (is_presel_event)
    {
      N_presel_events++;
      N_presel_stubs += nTotStubs;
      N_presel_stubs_0 += nStubs_0;
      N_presel_stubs_1 += nStubs_1;

      h_nTotStubs_presel->Fill(nTotStubs);
      h_nStubs_0_presel->Fill(nStubs_0);
      h_nStubs_1_presel->Fill(nStubs_1);
      h2_nStubs_presel->Fill(nStubs_0, nStubs_1);


      for (int k = 0; k < NMODULES; k++)
      {
        h_nstubsPerModule_presel->Fill(k, nstubs[k]);
      }
    }

    // S0 trackable and increasing requests in S1
    if (is_trackable_0 && fired_modules_1 >= 1)
    {
      N_trackable_0_onehit_1++;
      N_trackable_0_onehit_1_stubs += nTotStubs;
      N_trackable_0_onehit_1_stubs_0 += nStubs_0;
      N_trackable_0_onehit_1_stubs_1 += nStubs_1;

      h_nTotStubs_S0_trackable_S1_onehit->Fill(nTotStubs);
      h_nStubs_0_S0_trackable_S1_onehit->Fill(nStubs_0);
      h_nStubs_1_S0_trackable_S1_onehit->Fill(nStubs_1);
      h2_nStubs_S0_trackable_S1_onehit->Fill(nStubs_0, nStubs_1);

      for (int k = 0; k < NMODULES; k++)
      {
        h_nstubsPerModule_trackable_0_onehit_1->Fill(k, nstubs[k]);
      }
    }
    if (is_trackable_0 && fired_modules_1 >= 2)
    {
      N_trackable_0_twohit_1++;
    }
    if (is_trackable_0 && fired_modules_1 >= 3)
    {
      N_trackable_0_threehit_1++;
    }
    if (is_trackable_0 && fired_modules_1 >= 4)
    {
      N_trackable_0_fourhit_1++;
    }

    // S1 trackable and at least one hit in S0
    if (is_trackable_1 && fired_modules_0 >= 1)
    {
      N_trackable_1_onehit_0++;
      N_trackable_1_onehit_0_stubs += nTotStubs;
      N_trackable_1_onehit_0_stubs_0 += nStubs_0;
      N_trackable_1_onehit_0_stubs_1 += nStubs_1;

      h_nTotStubs_S1_trackable_S0_onehit->Fill(nTotStubs);
      h_nStubs_0_S1_trackable_S0_onehit->Fill(nStubs_0);
      h_nStubs_1_S1_trackable_S0_onehit->Fill(nStubs_1);
      h2_nStubs_S1_trackable_S0_onehit->Fill(nStubs_0, nStubs_1);

      for (int k = 0; k < NMODULES; k++)
      {
        h_nstubsPerModule_trackable_1_onehit_0->Fill(k, nstubs[k]);
      }
    }

    // at least a muon trackable in both stations
    bool is_trackable_event = is_trackable_0 && is_trackable_1;
    if (is_trackable_event)
    {
      N_trackable_events++;
      h_nTotStubs_trackable->Fill(nTotStubs);

      h_nStubs_0_trackable->Fill(nStubs_0);
      h_nStubs_1_trackable->Fill(nStubs_1);
      h2_nStubs_trackable->Fill(nStubs_0, nStubs_1);

      for (int k = 0; k < NMODULES; k++)
      {
        h_nstubsPerModule_trackable->Fill(k, nstubs[k]);
      }
    }

    bool is_noise_event = is_noise_0 && is_noise_1;
    if (is_noise_event)
    {
      N_noise_events++;
      h_nTotStubs_noise->Fill(nTotStubs);
    }

    bool is_missing_0_event = is_noise_0 && is_trackable_1;
    if (is_missing_0_event)
    {
      N_missing_0_events++;
      h_nTotStubs_missing_0->Fill(nTotStubs);
    }

    bool is_missing_1_event = is_trackable_0 && is_noise_1;
    if (is_missing_1_event)
    {
      N_missing_1_events++;
      h_nTotStubs_missing_1->Fill(nTotStubs);
    }

    bool is_not_trackable_event = !is_trackable_0 || !is_trackable_1;
    if (is_not_trackable_event)
    {
      N_not_trackable_events++;
      h_nTotStubs_not_trackable->Fill(nTotStubs);
    }

    bool is_presel_not_trackable_event = is_presel_event && is_not_trackable_event;
    if (is_presel_not_trackable_event)
    {
      N_presel_not_trackable_events++;
      h_nTotStubs_presel_not_trackable->Fill(nTotStubs);
    }

    bool is_fired12_event = is_fired6_0 && is_fired6_1;
    if (is_fired12_event)
    {
      N_fired12_events++;
      h_nTotStubs_fired12->Fill(nTotStubs);
      h_nStubs_0_fired12->Fill(nStubs_0);
      h_nStubs_1_fired12->Fill(nStubs_1);
      h2_nStubs_fired12->Fill(nStubs_0, nStubs_1);

      for (int k = 0; k < NMODULES; k++)
      {
        h_nstubsPerModule_fired12->Fill(k, nstubs[k]);
      }
    }

    bool is_passing_golden_event = is_passing_golden_0 && is_passing_golden_1;
    if (is_passing_golden_event)
    {

      N_passing_1mu_golden_events++;
    }

    //////////////////////////////3 stations

    bool is_passing_golden_2 = nStubs_2 == 6 && S_2.fired_xy_modules == 4 && S_2.fired_uv_modules == 2;
    if (is_passing_golden_2)
    {
      N_passing_1mu_golden_2++;
    }

    bool is_passing_golden_event_1 = is_passing_golden_1 && is_passing_golden_2;
    if (is_passing_golden_event_1)
    {

      N_passing_1mu_golden_events_1++;
    }

    bool is_presel_passing_golden_S0_event = is_presel_event && is_passing_golden_0;
    if (is_presel_passing_golden_S0_event)
    {
      N_presel_passing_golden_S0_events++;

      // plots LUMI
      // to test the absolute normalisation
      //
      // loop on all stubs in the two stations
      double X(-1);
      double Y(-1);
      bool foundX = false;
      bool foundY = false;
      for (int j = 0; j < nTotStubs; ++j)
      {
        // imod=(0:5) in station 0; imod=(6:11) in station 1
        int imod = Link->at(j);
        if (imod == 4)
        {
          X = LocalX->at(j);
          foundX = true;
          continue;
        }
        if (imod == 5)
        {
          Y = LocalX->at(j);
          foundY = true;
          continue;
        }
        if (foundX && foundY)
          break;
      }

      h2_LocalXY_presel_passing_golden->Fill(X, Y);

      if (iskim == 101)
        ok_skim = true;
    }

    bool is_presel_passing_golden_S1_event = is_presel_event && is_passing_golden_1;
    if (is_presel_passing_golden_S1_event)
    {
      N_presel_passing_golden_S1_events++;
    }

    //*************** efficiency 2 mu ***********************
    bool is_presel_passing_2mu_golden_S0_event = is_presel_event && is_passing_2mu_golden_0;
    if (is_presel_passing_2mu_golden_S0_event)
    {
      N_presel_passing_2mu_golden_S0_events++;
    }
    bool is_presel_passing_2mu_golden_S1_event = is_presel_event && is_passing_2mu_golden_1;
    if (is_presel_passing_2mu_golden_S1_event)
    {
      N_presel_passing_2mu_golden_S1_events++;
    }
    bool is_presel_passing_2mu_golden_S0_trackable_S1_event = is_presel_passing_2mu_golden_S0_event && is_passing_2mu_trackable_1;
    if (is_presel_passing_2mu_golden_S0_trackable_S1_event)
    {
      N_presel_passing_2mu_golden_S0_trackable_S1_events++;
    }
    bool is_presel_passing_2mu_golden_S1_trackable_S0_event = is_presel_passing_2mu_golden_S1_event && is_passing_2mu_trackable_0;
    if (is_presel_passing_2mu_golden_S1_trackable_S0_event)
    {
      N_presel_passing_2mu_golden_S1_trackable_S0_events++;
    }

    bool is_passing_2mu_golden_event = is_passing_2mu_golden_0 && is_passing_2mu_golden_1;
    if (is_passing_2mu_golden_event)
    {

      N_passing_2mu_golden_events++;
    }

    bool is_passing_2mu_golden_event_1 = is_passing_2mu_golden_1 && is_passing_2mu_golden_2;
    if (is_passing_2mu_golden_event_1)
    {

      N_passing_2mu_golden_events_1++;
    }

    //*************** efficiency 3 mu ***********************
    bool is_presel_passing_3mu_golden_S0_event = is_presel_event && is_passing_3mu_golden_0;
    if (is_presel_passing_3mu_golden_S0_event)
    {
      N_presel_passing_3mu_golden_S0_events++;
    }
    bool is_presel_passing_3mu_golden_S1_event = is_presel_event && is_passing_3mu_golden_1;
    if (is_presel_passing_3mu_golden_S1_event)
    {
      N_presel_passing_3mu_golden_S1_events++;
    }
    bool is_presel_passing_3mu_golden_S0_trackable_S1_event = is_presel_passing_3mu_golden_S0_event && is_passing_3mu_trackable_1;
    if (is_presel_passing_3mu_golden_S0_trackable_S1_event)
    {
      N_presel_passing_3mu_golden_S0_trackable_S1_events++;
    }
    bool is_presel_passing_3mu_golden_S1_trackable_S0_event = is_presel_passing_3mu_golden_S1_event && is_passing_3mu_trackable_0;
    if (is_presel_passing_3mu_golden_S1_trackable_S0_event)
    {
      N_presel_passing_3mu_golden_S1_trackable_S0_events++;
    }

    bool is_passing_3mu_golden_event = is_passing_3mu_golden_0 && is_passing_3mu_golden_1;
    if (is_passing_3mu_golden_event)
    {

      N_passing_3mu_golden_events++;
    }

    bool is_passing_3mu_golden_event_1 = is_passing_3mu_golden_1 && is_passing_3mu_golden_2;
    if (is_passing_3mu_golden_event_1)
    {

      N_passing_3mu_golden_events_1++;
    }

    //*************** efficiency 4 mu ***********************
    bool is_presel_passing_4mu_golden_S0_event = is_presel_event && is_passing_4mu_golden_0;
    if (is_presel_passing_4mu_golden_S0_event)
    {
      N_presel_passing_4mu_golden_S0_events++;
    }
    bool is_presel_passing_4mu_golden_S1_event = is_presel_event && is_passing_4mu_golden_1;
    if (is_presel_passing_4mu_golden_S1_event)
    {
      N_presel_passing_4mu_golden_S1_events++;
    }
    bool is_presel_passing_4mu_golden_S0_trackable_S1_event = is_presel_passing_4mu_golden_S0_event && is_passing_4mu_trackable_1;
    if (is_presel_passing_4mu_golden_S0_trackable_S1_event)
    {
      N_presel_passing_4mu_golden_S0_trackable_S1_events++;
    }
    bool is_presel_passing_4mu_golden_S1_trackable_S0_event = is_presel_passing_4mu_golden_S1_event && is_passing_4mu_trackable_0;
    if (is_presel_passing_4mu_golden_S1_trackable_S0_event)
    {
      N_presel_passing_4mu_golden_S1_trackable_S0_events++;
    }

    bool is_passing_4mu_golden_event = is_passing_4mu_golden_0 && is_passing_4mu_golden_1;
    if (is_passing_4mu_golden_event)
    {

      N_passing_4mu_golden_events++;
    }
    bool is_passing_4mu_golden_event_1 = is_passing_4mu_golden_1 && is_passing_4mu_golden_2;
    if (is_passing_4mu_golden_event_1)
    {

      N_passing_4mu_golden_events_1++;
    }
    //***************

    bool is_passing_golden_S0_zero_S1_event = is_passing_golden_0 && is_zero_1;
    if (is_passing_golden_S0_zero_S1_event)
      N_passing_golden_S0_zero_S1_events++;
    bool is_passing_golden_S1_zero_S0_event = is_passing_golden_1 && is_zero_0;
    if (is_passing_golden_S1_zero_S0_event)
      N_passing_golden_S1_zero_S0_events++;

    bool is_passing_golden_S0_trackable_S1_event = is_passing_golden_0 && is_trackable_1;
    if (is_passing_golden_S0_trackable_S1_event)
    {
      N_passing_golden_S0_trackable_S1_events++;
    }
    bool is_passing_golden_S1_trackable_S0_event = is_passing_golden_1 && is_trackable_0;
    if (is_passing_golden_S1_trackable_S0_event)
    {
      N_passing_golden_S1_trackable_S0_events++;
    }

    bool is_passing_1mu_event = is_single_clean_0 && is_passing_clean_1;
    if (is_passing_1mu_event)
    {
      N_passing_1mu_events++;
      h_nTotStubs_passing_1mu->Fill(nTotStubs);
    }
    bool is_passing_2mu_event = is_pileup_2mu_0 && is_pileup_2mu_1;
    if (is_passing_2mu_event)
    {
      N_passing_2mu_events++;
      h_nTotStubs_passing_2mu->Fill(nTotStubs);
    }
    bool is_passing_3mu_event = is_pileup_3mu_0 && is_pileup_3mu_1;
    if (is_passing_3mu_event)
    {
      N_passing_3mu_events++;
      h_nTotStubs_passing_3mu->Fill(nTotStubs);
    }
    bool is_passing_4mu_event = is_pileup_4mu_0 && is_pileup_4mu_1;
    if (is_passing_4mu_event)
    {
      N_passing_4mu_events++;
      h_nTotStubs_passing_4mu->Fill(nTotStubs);
    }

    bool is_single_cand_event = is_single_cand_0 && is_2tracks_1;
    if (is_single_cand_event)
    {
      N_single_cand_events++;
      h_nTotStubs_single_cand->Fill(nTotStubs);
    }
    /////////////////////// first pair: stations 0-1///////
    bool is_single_clean_event_0 = is_single_clean_0 && is_2tracks_1 && S_1.multifired_last_modules >= 1 && S_1.multifired_first_modules >= 1; // additional condition

    if (is_single_clean_event_0)
    {  
      N_single_clean_events++;
      h_nTotStubs_single_clean->Fill(nTotStubs);
      OfflineSkimMode.StationIndex = 0; // Station 0
      OfflineSkimMode.SingleMu = true;
      if (iskim == 1)
        ok_skim = true;
      if (iskim == 1001)
        ok_skim_MC = true;

      N_single_clean_stubs += nTotStubs;
      N_single_clean_stubs_0 += nStubs_0;
      N_single_clean_stubs_1 += nStubs_1;

      h_nStubs_0_single_clean->Fill(nStubs_0);
      h_nStubs_1_single_clean->Fill(nStubs_1);
      h2_nStubs_single_clean->Fill(nStubs_0, nStubs_1);


      for (int k = 0; k < NMODULES; k++)
      {
        h_nstubsPerModule_single_clean->Fill(k, nstubs[k]);
      }
    }

    //////////////////////////////////
    ///////////////////////////////////7
    /////second station///////////////

    bool is_single_clean_event_1 = is_single_clean_1 && is_2tracks_2 && S_2.multifired_last_modules >= 1 && S_2.multifired_first_modules >= 1; // additional condition
    if (is_single_clean_event_1)
    {

      N_single_clean_events_1++;
      h_nTotStubs_single_clean_2->Fill(nTotStubs);
      OfflineSkimMode.StationIndex = 1; // Station 0
      OfflineSkimMode.SingleMu = true;
      if (iskim == 1)
        ok_skim = true;
    }

    bool is_signal_plus_pileup_2 =
        S_2.multifired_xy_modules == 4 &&
        ((S_2.threefired_xy_modules >= 2 && S_2.nhits_uv >= 3) ||
         (S_2.threefired_xy_modules >= 3 && S_2.nhits_uv >= 2));
    if (is_signal_plus_pileup_2)
      N_signal_plus_pileup_2++;

    bool is_pileup_2mu_event_1 = is_pileup_2mu_1 && is_signal_plus_pileup_2 && S_2.threefired_last_modules >= 1 && S_2.threefired_first_modules >= 1;
    if (is_pileup_2mu_event_1)
    {
      N_pileup_2mu_events_1++;
    }

    bool is_signal_plus_tripileup_2 =
        S_2.threefired_xy_modules == 4 &&
        ((S_2.fourfired_xy_modules >= 2 && S_2.nhits_uv >= 4) ||
         (S_2.fourfired_xy_modules >= 3 && S_2.nhits_uv >= 3));
    if (is_signal_plus_tripileup_2)
      N_signal_plus_tripileup_2++;

    bool is_pileup_3mu_event_1 = is_pileup_3mu_1 && is_signal_plus_tripileup_2 && S_2.fourfired_last_modules >= 1 && S_2.fourfired_first_modules >= 1;
    if (is_pileup_3mu_event_1)
    {
      N_pileup_3mu_events_1++;
      // h_nTotStubs_pileup_3mu->Fill(nTotStubs);
    }

    bool is_signal_plus_fourpileup_2 =
        S_2.fourfired_xy_modules == 4 &&
        ((S_2.fivefired_xy_modules >= 2 && S_2.nhits_uv >= 5) ||
         (S_2.fivefired_xy_modules >= 3 && S_2.nhits_uv >= 4));
    if (is_signal_plus_fourpileup_2)
      N_signal_plus_fourpileup_2++;

    bool is_pileup_4mu_event_1 = is_pileup_4mu_1 && is_signal_plus_fourpileup_2 && S_2.fivefired_last_modules >= 1 && S_2.fivefired_first_modules >= 1;
    if (is_pileup_4mu_event_1)
    {
      N_pileup_4mu_events_1++;
      // h_nTotStubs_pileup_4mu->Fill(nTotStubs);
    }

    bool is_pileup234_skim_event_1 = is_pileup_2mu_event_1 || is_pileup_3mu_event_1 || is_pileup_4mu_event_1;
    if (is_pileup234_skim_event_1)
    {
      N_pileup234_skim_events_1++;
      // h_nTotStubs_pileup234_skim->Fill(nTotStubs);
      if (iskim == 60)
        ok_skim = true;

      // N_pileup234_skim_stubs += nTotStubs;
      // N_pileup234_skim_stubs_0 += nStubs_0;
      // N_pileup234_skim_stubs_1 += nStubs_1;

      // h_nStubs_0_pileup234_skim->Fill(nStubs_0);
      // h_nStubs_1_pileup234_skim->Fill(nStubs_1);
      // h2_nStubs_pileup234_skim->Fill(nStubs_0, nStubs_1);

  

      // for (int k = 0; k < NMODULES; k++)
      // {
      //   h_nstubsPerModule_pileup234_skim->Fill(k, nstubs[k]);
      // }
    }

    bool is_fired12_plus_single_clean_event = is_fired6_plus_single_clean_0 && is_fired6_2tracks_1;
    if (is_fired12_plus_single_clean_event)
    {
      N_fired12_plus_single_clean_events++;
      h_nTotStubs_fired12_plus_single_clean->Fill(nTotStubs);

      h_nStubs_0_fired12_plus_single_clean->Fill(nStubs_0);
      h_nStubs_1_fired12_plus_single_clean->Fill(nStubs_1);
      h2_nStubs_fired12_plus_single_clean->Fill(nStubs_0, nStubs_1);

      for (int k = 0; k < NMODULES; k++)
      {
        h_nstubsPerModule_fired12_plus_single_clean->Fill(k, nstubs[k]);
      }
    }

    bool is_umsel_event = nStubs_0 >= 5 && nStubs_1 >= 5 && (nStubs_1 - nStubs_0) >= 5;
    if (is_umsel_event)
    {
      N_umsel_events++;
      h_nTotStubs_umsel->Fill(nTotStubs);
    }

    bool is_golden_event = is_passing_golden_0 && is_golden_1;
    if (is_golden_event)
    {
      N_golden_events++;
      h_nTotStubs_golden->Fill(nTotStubs);
    }

    bool is_pileup_any_event = is_pileup_any_0 && is_signal_plus_pileup_1;
    if (is_pileup_any_event)
    {
      N_pileup_any_events++;
      h_nTotStubs_pileup_any->Fill(nTotStubs);
    }
    bool is_pileup_2mu_event = is_pileup_2mu_0 && is_signal_plus_pileup_1 && S_1.threefired_last_modules >= 1 && S_1.threefired_first_modules >= 1;
    if (is_pileup_2mu_event)
    {
      N_pileup_2mu_events++;
      h_nTotStubs_pileup_2mu->Fill(nTotStubs);
    }
    bool is_pileup_3mu_event = is_pileup_3mu_0 && is_signal_plus_tripileup_1 && S_1.fourfired_last_modules >= 1 && S_1.fourfired_first_modules >= 1;
    if (is_pileup_3mu_event)
    {
      N_pileup_3mu_events++;
      h_nTotStubs_pileup_3mu->Fill(nTotStubs);
    }
    bool is_pileup_4mu_event = is_pileup_4mu_0 && is_signal_plus_fourpileup_1 && S_1.fivefired_last_modules >= 1 && S_1.fivefired_first_modules >= 1;
    if (is_pileup_4mu_event)
    {
      N_pileup_4mu_events++;
      h_nTotStubs_pileup_4mu->Fill(nTotStubs);
    }
    bool is_pileup_many_event = is_pileup_many_0 && is_signal_plus_fourpileup_1;
    if (is_pileup_many_event)
    {
      N_pileup_many_events++;
      h_nTotStubs_pileup_many->Fill(nTotStubs);
    }

    bool is_pileup_skim_event = is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_many_event;
    if (is_pileup_skim_event)
    {
      N_pileup_skim_events++;
      h_nTotStubs_pileup_skim->Fill(nTotStubs);
      if (iskim == 40)
        ok_skim = true;

      N_pileup_skim_stubs += nTotStubs;
      N_pileup_skim_stubs_0 += nStubs_0;
      N_pileup_skim_stubs_1 += nStubs_1;

      h_nStubs_0_pileup_skim->Fill(nStubs_0);
      h_nStubs_1_pileup_skim->Fill(nStubs_1);
      h2_nStubs_pileup_skim->Fill(nStubs_0, nStubs_1);

    }

    bool is_pileup234_skim_event = is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_4mu_event;
    if (is_pileup234_skim_event)
    {
      N_pileup234_skim_events_1++;
      h_nTotStubs_pileup234_skim->Fill(nTotStubs);
      if (iskim == 60)
        ok_skim = true;

      N_pileup234_skim_stubs += nTotStubs;
      N_pileup234_skim_stubs_0 += nStubs_0;
      N_pileup234_skim_stubs_1 += nStubs_1;

      h_nStubs_0_pileup234_skim->Fill(nStubs_0);
      h_nStubs_1_pileup234_skim->Fill(nStubs_1);
      h2_nStubs_pileup234_skim->Fill(nStubs_0, nStubs_1);


      for (int k = 0; k < NMODULES; k++)
      {
        h_nstubsPerModule_pileup234_skim->Fill(k, nstubs[k]);
      }
    }

    bool is_loose_event = is_single_clean_event_0 || is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_many_event;
    if (is_loose_event)
    {
      N_loose_events++;
      if (iskim == 41)
        ok_skim = true;
    }

    bool is_loose234_event = is_single_clean_event_0 || is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_4mu_event;
    if (is_loose234_event)
    {
      N_loose234_events++;
      if (iskim == 61)
        ok_skim = true;
    }

    if (ok_skim)
    {
      o_tree->Fill();
      o_nevt++;

      N_stubs_skim += nTotStubs;
    }
    else if (ok_skim_MC)
    {
      o_tree->Fill();
      o_nevt++;

      N_stubs_skim += nTotStubs;
    }

  } //  for(Long64_t i = 0; i < n_entries; i++)

  cout << " N events: " << n_entries << endl;
  cout << " N empty events: " << N_empty << endl;
  cout << " N diff size events: " << N_differentSize << endl;

  cout << " N_single_clean_events: " << N_single_clean_events << endl;
  cout << " N_single_clean_events_1: " << N_single_clean_events_1 << endl;

  printEff("offline+online single_muon_interaction_0",
           N_single_muon_offline_online_0,
           N_single_clean_events);

  printEff("offline+online single_muon_interaction_1",
           N_single_muon_offline_online_1,
           N_single_clean_events_1);

  printEff("offline+online pileup_muon_interaction_0",
           N_pileup_muon_offline_online_0,
           N_pileup_2mu_events + N_pileup_3mu_events + N_pileup_4mu_events);

  printEff("offline+online pileup_muon_interaction_1",
           N_pileup_muon_offline_online_1,
           N_pileup_2mu_events_1 + N_pileup_3mu_events_1 + N_pileup_4mu_events_1);

  printEff("offline+online single_passing_muon_0",
           N_single_passing_muon_offline_online_0,
           N_passing_1mu_golden_events);

  printEff("offline+online single_passing_muon_1",
           N_single_passing_muon_offline_online_1,
           N_passing_1mu_golden_events_1);

  printEff("offline+online pileup_passing_muon_0",
           N_pileup_passing_muon_offline_online_0,
           N_passing_2mu_golden_0 + N_passing_3mu_golden_0 + N_passing_4mu_golden_0);

  printEff("offline+online pileup_passing_muon_1",
           N_pileup_passing_muon_offline_online_1,
           N_passing_2mu_golden_events_1 + N_passing_3mu_golden_events_1 + N_passing_4mu_golden_events_1);



  cout << "\n offline online single_muon_interaction_0:" << N_single_muon_offline_online_0 << endl; // N_single_muon_offline_online /
  cout << "\n offline st0:" << N_single_clean_events << endl;
  cout << "\n eff:" << double(N_single_muon_offline_online_0) / double(N_single_clean_events) << endl;

  cout << "--------------------------------------" << endl;

  cout << "\n offline online single_muon_interaction_1:" << N_single_muon_offline_online_1 << endl; // N_single_muon_offline_online /
  cout << "\n offline st1:" << N_single_clean_events_1 << endl;
  cout << "\n eff:" << double(N_single_muon_offline_online_1) / double(N_single_clean_events_1) << endl;

  cout << "--------------------------------------" << endl;

  cout << "\n offline online pileup_muon_interaction_0:" << N_pileup_muon_offline_online_0 << endl; // N_single_muon_offline_online /
  cout << "\n offline st0:" << N_pileup_2mu_events + N_pileup_3mu_events + N_pileup_4mu_events << endl;
  cout << "\n eff:" << double(N_pileup_muon_offline_online_0) / double(N_pileup_2mu_events + N_pileup_3mu_events + N_pileup_4mu_events) << endl;

  cout << "--------------------------------------" << endl;

  cout << "\n offline online pileup_muon_interaction_1:" << N_pileup_muon_offline_online_1 << endl; // N_single_muon_offline_online /
  cout << "\n offline st1:" << N_pileup_2mu_events_1 + N_pileup_3mu_events_1 + N_pileup_4mu_events_1 << endl;
  cout << "\n eff:" << double(N_pileup_muon_offline_online_1) / double(N_pileup_2mu_events_1 + N_pileup_3mu_events_1 + N_pileup_4mu_events_1) << endl;

  cout << "--------------------------------------" << endl;

  cout << "\n offline online single_passing_muon_0:" << N_single_passing_muon_offline_online_0 << endl; // N_single_muon_offline_online /
  cout << "\n offline st0:" << N_passing_1mu_golden_events << endl;
  cout << "\n eff:" << double(N_single_passing_muon_offline_online_0) / double(N_passing_1mu_golden_events) << endl;

  cout << "--------------------------------------" << endl;

  cout << "\n offline online single_passing_muon_1 (GOLDEN):" << N_single_passing_muon_offline_online_1 << endl; // N_single_muon_offline_online /
  cout << "\n offline st1" << N_passing_1mu_golden_events_1 << endl;
  cout << "\n eff:" << double(N_single_passing_muon_offline_online_1) / double(N_passing_1mu_golden_events_1) << endl;

  cout << "--------------------------------------" << endl;

  cout << "\n offline online pileup_passing_muon_0:" << N_pileup_passing_muon_offline_online_0 << endl; // N_single_muon_offline_online /
  cout << "\n offline st0:" << N_passing_2mu_golden_0 + N_passing_3mu_golden_0 + N_passing_4mu_golden_0 << endl;
  cout << "\n eff:" << double(N_pileup_passing_muon_offline_online_0) / double(N_passing_2mu_golden_0 + N_passing_3mu_golden_0 + N_passing_4mu_golden_0) << endl;

  cout << "--------------------------------------" << endl;

  cout << "\n offline online pileup_passing_muon_1 (GOLDEN):" << N_pileup_passing_muon_offline_online_1 << endl; // N_single_muon_offline_online /
  cout << "\n offline st1:" << N_passing_2mu_golden_events_1 + N_passing_3mu_golden_events_1 + N_passing_4mu_golden_events_1 << endl;
  cout << "\n eff:" << double(N_pileup_passing_muon_offline_online_1) / double(N_passing_2mu_golden_events_1 + N_passing_3mu_golden_events_1 + N_passing_4mu_golden_events_1) << endl;

  cout << "\n"
       << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\n Station 0:" << endl;
  cout << "Number of events without any hit in S_0 = " << endl;
  cout << "\t" << N_zero_0 << endl;
  cout << "Number of Empty or Noise events in S_0 = " << endl;
  cout << "\t" << N_noise_0 << endl;
  cout << "Number of events with >=1 fired S_0 module = " << endl;
  cout << "\t" << N_onehit_0 << endl;
  cout << "Number of events with >=2 fired S_0 modules = " << endl;
  cout << "\t" << N_twohit_0 << endl;
  cout << "Number of events with >=3 fired S_0 modules = " << endl;
  cout << "\t" << N_threehit_0 << endl;
  cout << "Number of events with >=4 fired S_0 modules = " << endl;
  cout << "\t" << N_fourhit_0 << endl;
  cout << "Number of events with all 6 fired S_0 modules = " << endl;
  cout << "\t" << N_fired6_0 << endl;
  cout << "Number of trackable events (all 4 XY modules and at least 1 UV) = " << endl;
  cout << "\t " << N_trackable_0 << endl;
  cout << "Number of candidate single mu events in S0 = " << endl;
  cout << "\t " << N_single_cand_0 << endl;
  cout << "Number of clean single mu events in S_0 = " << endl;
  cout << "\t " << N_single_clean_0 << endl;
  cout << "Number of golden 1mu passing events in Station 0 (all 6 modules with one hit) = " << endl;
  cout << "\t" << N_passing_1mu_golden_0 << endl;
  cout << "Number of golden 2mu passing events in Station 0 (all 6 modules with two hits) = " << endl;
  cout << "\t" << N_passing_2mu_golden_0 << endl;
  cout << "Number of golden 3mu passing events in Station 0 (all 6 modules with three hits) = " << endl;
  cout << "\t" << N_passing_3mu_golden_0 << endl;
  cout << "Number of golden 4mu passing events in Station 0 (all 6 modules with four hits) = " << endl;
  cout << "\t" << N_passing_4mu_golden_0 << endl;
  cout << "Number of (any) pileup events in S_0 = " << endl;
  cout << "\t " << N_pileup_any_0 << endl;
  cout << "Number of 2mu pileup events in S_0 = " << endl;
  cout << "\t " << N_pileup_2mu_0 << endl;
  cout << "Number of 3mu pileup events in S_0 = " << endl;
  cout << "\t " << N_pileup_3mu_0 << endl;
  cout << "Number of 4mu pileup events in S_0 = " << endl;
  cout << "\t " << N_pileup_4mu_0 << endl;
  cout << "Number of >=4mu pileup events in S_0 = " << endl;
  cout << "\t " << N_pileup_many_0 << endl;

  cout << "\n Station 1:" << endl;
  cout << "Number of events without any hit in S_1 = " << endl;
  cout << "\t" << N_zero_1 << endl;
  cout << "Number of Empty or Noise events in S_1 = " << endl;
  cout << "\t" << N_noise_1 << endl;
  cout << "Number of events with all 6 fired S_1 modules = " << endl;
  cout << "\t" << N_fired6_1 << endl;
  cout << "Number of golden 1mu passing events in Station 1 (all 6 modules with one hit) = " << endl;
  cout << "\t" << N_passing_1mu_golden_1 << endl;
  cout << "Number of golden 2mu passing events in Station 1 (all 6 modules with two hits) = " << endl;
  cout << "\t" << N_passing_2mu_golden_1 << endl;
  cout << "Number of golden 3mu passing events in Station 1 (all 6 modules with three hits) = " << endl;
  cout << "\t" << N_passing_3mu_golden_1 << endl;
  cout << "Number of golden 4mu passing events in Station 1 (all 6 modules with four hits) = " << endl;
  cout << "\t" << N_passing_4mu_golden_1 << endl;
  cout << "Number of candidate single mu events in S_1 = " << endl;
  cout << "\t " << N_passing_cand_1 << endl;
  cout << "Number of clean single mu events in S_1 = " << endl;
  cout << "\t " << N_passing_clean_1 << endl;
  cout << "Number of interaction candidates (2 tracks) in S_1 = " << endl;
  cout << "\t " << N_2tracks_1 << endl;
  cout << "Number of golden events in S_1 (2 hits in all modules in S_1)= " << endl;
  cout << "\t" << N_golden_1 << endl;
  cout << "Number of signal+pileup candidates in S_1 = " << endl;
  cout << "\t " << N_signal_plus_pileup_1 << endl;

  cout << "\n EVENTS:   " << endl;
  cout << "Number all input events = " << endl;
  cout << "\t" << N_input_events << endl;
  cout << "Number of events without any hit in both S_0 and S_1 = " << endl;
  cout << "\t" << N_zero_events << endl;
  cout << "Number of Noise events (empty or <=1 module in both stations) = " << endl;
  cout << "\t " << N_noise_events << endl;
  cout << "Number of events missing S_0 (empty or <=1 module) and trackable in S_1 = " << endl;
  cout << "\t " << N_missing_0_events << endl;
  cout << "Number of events missing S_1 (empty or <=1 module) and trackable in S_0 = " << endl;
  cout << "\t " << N_missing_1_events << endl;
  cout << "Number of Preselected events = " << endl;
  cout << "\t " << N_presel_events << endl;
  //
  cout << "Number of trackable S_0 events with >=1 fired S_1 module = " << endl;
  cout << "\t" << N_trackable_0_onehit_1 << endl;
  cout << "Number of trackable S_0 events with >=2 fired S_1 module = " << endl;
  cout << "\t" << N_trackable_0_twohit_1 << endl;
  cout << "Number of trackable S_0 events with >=3 fired S_1 module = " << endl;
  cout << "\t" << N_trackable_0_threehit_1 << endl;
  cout << "Number of trackable S_0 events with >=4 fired S_1 module = " << endl;
  cout << "\t" << N_trackable_0_fourhit_1 << endl;
  cout << "Number of trackable events = " << endl; // both stations trackable
  cout << "\t " << N_trackable_events << endl;

  cout << "Number of events with ALL 12 fired modules = " << endl;
  cout << "\t" << N_fired12_events << endl;

  cout << "Number of clean 1mu passing events = " << endl;
  cout << "\t " << N_passing_1mu_events << endl;

  cout << "Number of passing 1mu golden S_0 in preselected events = " << endl;
  cout << "\t " << N_presel_passing_golden_S0_events << endl;
  cout << "Number of passing 1mu golden S_0 and trackable in S_1 = " << endl;
  cout << "\t " << N_passing_golden_S0_trackable_S1_events << endl;
  cout << "\tNumber of passing 1mu golden S_0 and zero hits in S_1 = " << endl;
  cout << "\t\t " << N_passing_golden_S0_zero_S1_events << endl;
  //
  cout << "Number of passing 1mu golden S_1 in preselected events = " << endl;
  cout << "\t " << N_presel_passing_golden_S1_events << endl;
  cout << "Number of passing 1mu golden S_1 and trackable in S_0 = " << endl;
  cout << "\t " << N_passing_golden_S1_trackable_S0_events << endl;
  cout << "\tNumber of passing 1mu golden S_1 and zero hits in S_0 = " << endl;
  cout << "\t\t " << N_passing_golden_S1_zero_S0_events << endl;

  cout << "Number of golden 1mu passing events (all 12 modules with one hit) = " << endl;
  cout << "\t" << N_passing_1mu_golden_events << endl;
  cout << "Number of golden 2mu passing events (all 12 modules with two hits) = " << endl;
  cout << "\t" << N_passing_2mu_golden_events << endl;
  cout << "Number of golden 3mu passing events (all 12 modules with three hits) = " << endl;
  cout << "\t" << N_passing_3mu_golden_events << endl;
  cout << "Number of golden 4mu passing events (all 12 modules with four hits) = " << endl;
  cout << "\t" << N_passing_4mu_golden_events << endl;

  cout << "Number of candidate Single mu interaction events = " << endl;
  cout << "\t " << N_single_cand_events << endl;
  cout << "Number of clean Single mu interaction events = " << endl;
  cout << "\t " << N_single_clean_events << endl;
  cout << "Number of golden Single mu interaction events = " << endl;
  cout << "\t " << N_golden_events << endl;
  cout << "Number of any pileup interaction events = " << endl;
  cout << "\t " << N_pileup_any_events << endl;
  cout << "Number of 2mu pileup interaction events = " << endl;
  cout << "\t " << N_pileup_2mu_events << endl;
  cout << "Number of 3mu pileup interaction events = " << endl;
  cout << "\t " << N_pileup_3mu_events << endl;
  cout << "Number of 4mu pileup interaction events = " << endl;
  cout << "\t " << N_pileup_4mu_events << endl;
  cout << "Number of >=4mu pileup interaction events = " << endl;
  cout << "\t " << N_pileup_many_events << endl;
  cout << "Total number of pileup (2,3,4) interaction events = " << endl;
  cout << "\t " << N_pileup234_skim_events << endl;
  cout << "Total number of pileup (2,3,>=4) interaction events = " << endl;
  cout << "\t " << N_pileup_skim_events << endl;
  cout << "Total number of Single + Pileup (2,3,4) interaction events = " << endl;
  cout << "\t " << N_loose234_events << endl;
  cout << "Total number of Single + Pileup (2,3,>=4) interaction events = " << endl;
  cout << "\t " << N_loose_events << endl;
  cout << "Number of events with 1mu (12stubs) overlapping Single Clean events = " << endl;
  cout << "\t" << N_fired12_plus_single_clean_events << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  cout << "N INPUT Events    = " << N_input_events << endl;
  cout << "N input stubs S_0 = " << N_input_stubs_0 << endl;
  cout << "N input stubs S_1 = " << N_input_stubs_1 << endl;
  cout << "N input stubs tot = " << N_input_stubs << endl;
  cout << endl;
  cout << "Fraction missing S_0             = " << double(N_missing_0_events) / double(N_input_events) << endl;
  cout << "Fraction missing S_1             = " << double(N_missing_1_events) / double(N_input_events) << endl;
  cout << "Fraction not trackable           = " << double(N_not_trackable_events) / double(N_input_events) << endl;

  cout << endl
       << "N PRESELECTED Events     = " << N_presel_events << endl;
  cout << "N stubs S_0 (presel.ev.) = " << N_presel_stubs_0 << endl;
  cout << "N stubs S_1 (presel.ev.) = " << N_presel_stubs_1 << endl;
  cout << "N stubs tot (presel.ev.) = " << N_presel_stubs << endl;

  cout << endl
       << "N PRESELECTED S_0 Trackable Events     = " << N_trackable_0_onehit_1 << endl;
  cout << "N stubs S_0 (presel.S0-trackable ev.) = " << N_trackable_0_onehit_1_stubs_0 << endl;
  cout << "N stubs S_1 (presel.S0-trackable ev.) = " << N_trackable_0_onehit_1_stubs_1 << endl;
  cout << "N stubs tot (presel.S0-trackable ev.) = " << N_trackable_0_onehit_1_stubs << endl;

  cout << endl
       << "N TRACKABLE Events       = " << N_trackable_events << endl;

  cout << "Fraction PRESELECTED / all       = " << double(N_presel_events) / double(N_input_events) << endl;
  cout << "Fraction PRESEL. S_0-Trackable / PRESEL = " << double(N_trackable_0_onehit_1) / double(N_presel_events) << endl;
  cout << "Fraction TRACKABLE / PRESEL      = " << double(N_trackable_events) / double(N_presel_events) << endl;

  // Average Station efficiencies
  pair<double, double> ineff_zero_S1 = binomial_eff(N_passing_golden_S0_zero_S1_events, N_passing_1mu_golden_0);
  pair<double, double> eff_S1 = binomial_eff(N_passing_golden_S0_trackable_S1_events, N_passing_1mu_golden_0);
  pair<double, double> eff_presel_S1 = binomial_eff(N_passing_golden_S0_trackable_S1_events, N_presel_passing_golden_S0_events);

  pair<double, double> ineff_zero_S0 = binomial_eff(N_passing_golden_S1_zero_S0_events, N_passing_1mu_golden_1);
  pair<double, double> eff_S0 = binomial_eff(N_passing_golden_S1_trackable_S0_events, N_passing_1mu_golden_1);
  pair<double, double> eff_presel_S0 = binomial_eff(N_passing_golden_S1_trackable_S0_events, N_presel_passing_golden_S1_events);

  cout << endl
       << "estimated Efficiency reconstruction S_1 (after presel) = " << eff_presel_S1.first << " +- " << eff_presel_S1.second << endl;
  cout << "estimated Efficiency reconstruction S_1 (all events)   = " << eff_S1.first << " +- " << eff_S1.second << endl;
  cout << "\t estimated inefficiency (zero hits in S_1)  = " << ineff_zero_S1.first << " +- " << ineff_zero_S1.second << endl;
  cout << "estimated Efficiency reconstruction S_0 (after presel) = " << eff_presel_S0.first << " +- " << eff_presel_S0.second << endl;
  cout << "estimated Efficiency reconstruction S_0 (all events)   = " << eff_S0.first << " +- " << eff_S0.second << endl;
  cout << "\t estimated inefficiency (zero hits in S_0)  = " << ineff_zero_S0.first << " +- " << ineff_zero_S0.second << endl;

  // plot station efficiency vs time
  h_effS1_vs_time->Divide(h_goldenS0_trackableS1_vs_time, h_goldenS0_vs_time, 1, 1, "B");
  h_effS0_vs_time->Divide(h_goldenS1_trackableS0_vs_time, h_goldenS1_vs_time, 1, 1, "B");

  h_presel_effS1_vs_time->Divide(h_goldenS0_trackableS1_vs_time, h_presel_goldenS0_vs_time, 1, 1, "B");
  h_presel_effS0_vs_time->Divide(h_goldenS1_trackableS0_vs_time, h_presel_goldenS1_vs_time, 1, 1, "B");

  // Average 2,3,4-mu station efficiencies (Trackable definition)
  pair<double, double> eff_presel_2mu_S1 = binomial_eff(N_presel_passing_2mu_golden_S0_trackable_S1_events, N_presel_passing_2mu_golden_S0_events);
  pair<double, double> eff_presel_2mu_S0 = binomial_eff(N_presel_passing_2mu_golden_S1_trackable_S0_events, N_presel_passing_2mu_golden_S1_events);
  //
  pair<double, double> eff_presel_3mu_S1 = binomial_eff(N_presel_passing_3mu_golden_S0_trackable_S1_events, N_presel_passing_3mu_golden_S0_events);
  pair<double, double> eff_presel_3mu_S0 = binomial_eff(N_presel_passing_3mu_golden_S1_trackable_S0_events, N_presel_passing_3mu_golden_S1_events);
  //
  pair<double, double> eff_presel_4mu_S1 = binomial_eff(N_presel_passing_4mu_golden_S0_trackable_S1_events, N_presel_passing_4mu_golden_S0_events);
  pair<double, double> eff_presel_4mu_S0 = binomial_eff(N_presel_passing_4mu_golden_S1_trackable_S0_events, N_presel_passing_4mu_golden_S1_events);
  // plot 2,3,4-mu station efficiencies  vs time (Trackable definition)
  h_presel_2mu_effS1_vs_time->Divide(h_presel_passing_2mu_goldenS0_trackableS1_vs_time, h_presel_passing_2mu_goldenS0_vs_time, 1, 1, "B");
  h_presel_2mu_effS0_vs_time->Divide(h_presel_passing_2mu_goldenS1_trackableS0_vs_time, h_presel_passing_2mu_goldenS1_vs_time, 1, 1, "B");
  h_presel_3mu_effS1_vs_time->Divide(h_presel_passing_3mu_goldenS0_trackableS1_vs_time, h_presel_passing_3mu_goldenS0_vs_time, 1, 1, "B");
  h_presel_3mu_effS0_vs_time->Divide(h_presel_passing_3mu_goldenS1_trackableS0_vs_time, h_presel_passing_3mu_goldenS1_vs_time, 1, 1, "B");
  h_presel_4mu_effS1_vs_time->Divide(h_presel_passing_4mu_goldenS0_trackableS1_vs_time, h_presel_passing_4mu_goldenS0_vs_time, 1, 1, "B");
  h_presel_4mu_effS0_vs_time->Divide(h_presel_passing_4mu_goldenS1_trackableS0_vs_time, h_presel_passing_4mu_goldenS1_vs_time, 1, 1, "B");

  // Average 2,3,4-mu station efficiencies (Golden definition)
  pair<double, double> eff_golden_1mu_S1 = binomial_eff(N_passing_1mu_golden_events, N_presel_passing_golden_S0_events);
  pair<double, double> eff_golden_1mu_S0 = binomial_eff(N_passing_1mu_golden_events, N_presel_passing_golden_S1_events);
  pair<double, double> eff_golden_2mu_S1 = binomial_eff(N_passing_2mu_golden_events, N_presel_passing_2mu_golden_S0_events);
  pair<double, double> eff_golden_2mu_S0 = binomial_eff(N_passing_2mu_golden_events, N_presel_passing_2mu_golden_S1_events);
  pair<double, double> eff_golden_3mu_S1 = binomial_eff(N_passing_3mu_golden_events, N_presel_passing_3mu_golden_S0_events);
  pair<double, double> eff_golden_3mu_S0 = binomial_eff(N_passing_3mu_golden_events, N_presel_passing_3mu_golden_S1_events);
  pair<double, double> eff_golden_4mu_S1 = binomial_eff(N_passing_4mu_golden_events, N_presel_passing_4mu_golden_S0_events);
  pair<double, double> eff_golden_4mu_S0 = binomial_eff(N_passing_4mu_golden_events, N_presel_passing_4mu_golden_S1_events);
  // plot 2,3,4-mu station efficiencies  vs time (Golden definition)
  h_golden_effS1_vs_time->Divide(h_passing_1mu_golden_vs_time, h_presel_goldenS0_vs_time, 1, 1, "B");
  h_golden_effS0_vs_time->Divide(h_passing_1mu_golden_vs_time, h_presel_goldenS1_vs_time, 1, 1, "B");
  h_golden_2mu_effS1_vs_time->Divide(h_passing_2mu_golden_vs_time, h_presel_passing_2mu_goldenS0_vs_time, 1, 1, "B");
  h_golden_2mu_effS0_vs_time->Divide(h_passing_2mu_golden_vs_time, h_presel_passing_2mu_goldenS1_vs_time, 1, 1, "B");
  h_golden_3mu_effS1_vs_time->Divide(h_passing_3mu_golden_vs_time, h_presel_passing_3mu_goldenS0_vs_time, 1, 1, "B");
  h_golden_3mu_effS0_vs_time->Divide(h_passing_3mu_golden_vs_time, h_presel_passing_3mu_goldenS1_vs_time, 1, 1, "B");
  h_golden_4mu_effS1_vs_time->Divide(h_passing_4mu_golden_vs_time, h_presel_passing_4mu_goldenS0_vs_time, 1, 1, "B");
  h_golden_4mu_effS0_vs_time->Divide(h_passing_4mu_golden_vs_time, h_presel_passing_4mu_goldenS1_vs_time, 1, 1, "B");
  ////////////////////////////////////////////////////
  ofstream ofeff1("relative_efficiency_S1_goodS0.txt");
  ofeff1 << "RELATIVE EFFICIENCY Station S_1 (in events with Golden Passing Mu in S_0)" << endl;
  ofeff1 << "#all_S0_passing_golden_events: " << N_passing_1mu_golden_0 << endl;
  ofeff1 << "#presel_S0_passing_golden_events: " << N_presel_passing_golden_S0_events << endl;
  ofeff1 << "#S0_passing_golden_&_S1_trackable: " << N_passing_golden_S0_trackable_S1_events << endl;
  ofeff1 << "#passing_golden_in_S0_&&_S1: " << N_passing_1mu_golden_events << endl;
  ofeff1 << "efficiency S1 trackable / all S0 passing golden = " << eff_S1.first << " +- " << eff_S1.second << endl;
  ofeff1 << "efficiency S1 trackable / presel S0 pass.golden = " << eff_presel_S1.first << " +- " << eff_presel_S1.second << endl;
  ofeff1 << "efficiency S1 golden    / presel S0 pass.golden = " << eff_golden_1mu_S1.first << " +- " << eff_golden_1mu_S1.second << endl;
  ofeff1 << "#S1_zero: " << N_passing_golden_S0_zero_S1_events << endl;
  ofeff1 << "fraction ZERO/all : " << ineff_zero_S1.first << " +- " << ineff_zero_S1.second << endl;
  //
  ofeff1 << "#presel_2mu_S0_passing_golden_events: " << N_presel_passing_2mu_golden_S0_events << endl;
  ofeff1 << "#2mu_S0_passing_golden_&_S1_trackable: " << N_presel_passing_2mu_golden_S0_trackable_S1_events << endl;
  ofeff1 << "#2mu_passing_golden_in_S0_&&_S1: " << N_passing_2mu_golden_events << endl;
  ofeff1 << "2mu efficiency S1 trackable / presel S0 passing golden = " << eff_presel_2mu_S1.first << " +- " << eff_presel_2mu_S1.second << endl;
  ofeff1 << "2mu efficiency S1 golden    / presel S0 pass.golden    = " << eff_golden_2mu_S1.first << " +- " << eff_golden_2mu_S1.second << endl;
  //
  ofeff1 << "#presel_3mu_S0_passing_golden_events: " << N_presel_passing_3mu_golden_S0_events << endl;
  ofeff1 << "#3mu_S0_passing_golden_&_S1_trackable: " << N_presel_passing_3mu_golden_S0_trackable_S1_events << endl;
  ofeff1 << "#3mu_passing_golden_in_S0_&&_S1: " << N_passing_3mu_golden_events << endl;
  ofeff1 << "3mu efficiency S1 trackable / presel S0 passing golden = " << eff_presel_3mu_S1.first << " +- " << eff_presel_3mu_S1.second << endl;
  ofeff1 << "3mu efficiency S1 golden    / presel S0 pass.golden    = " << eff_golden_3mu_S1.first << " +- " << eff_golden_3mu_S1.second << endl;
  //
  ofeff1 << "#presel_4mu_S0_passing_golden_events: " << N_presel_passing_4mu_golden_S0_events << endl;
  ofeff1 << "#4mu_S0_passing_golden_&_S1_trackable: " << N_presel_passing_4mu_golden_S0_trackable_S1_events << endl;
  ofeff1 << "#4mu_passing_golden_in_S0_&&_S1: " << N_passing_4mu_golden_events << endl;
  ofeff1 << "4mu efficiency S1 trackable / presel S0 passing golden = " << eff_presel_4mu_S1.first << " +- " << eff_presel_4mu_S1.second << endl;
  ofeff1 << "4mu efficiency S1 golden    / presel S0 pass.golden    = " << eff_golden_4mu_S1.first << " +- " << eff_golden_4mu_S1.second << endl;
  ofeff1.close();

  ofstream ofeff0("relative_efficiency_S0_goodS1.txt");
  ofeff0 << "RELATIVE EFFICIENCY Station S_0 (in events with Golden Passing Mu in S_1)" << endl;
  ofeff0 << "#all_S1_passing_golden_events: " << N_passing_1mu_golden_1 << endl;
  ofeff0 << "#presel_S1_passing_golden_events: " << N_presel_passing_golden_S1_events << endl;
  ofeff0 << "#S1_passing_golden_&_S0_trackable: " << N_passing_golden_S1_trackable_S0_events << endl;
  ofeff0 << "#passing_golden_in_S0_&&_S1: " << N_passing_1mu_golden_events << endl;
  ofeff0 << "efficiency S0 trackable / all S1 passing golden = " << eff_S0.first << " +- " << eff_S0.second << endl;
  ofeff0 << "efficiency S0 trackable / presel S1 pass.golden = " << eff_presel_S0.first << " +- " << eff_presel_S0.second << endl;
  ofeff0 << "efficiency S0 golden    / presel S1 pass.golden = " << eff_golden_1mu_S0.first << " +- " << eff_golden_1mu_S0.second << endl;
  ofeff0 << "#S0_zero: " << N_passing_golden_S1_zero_S0_events << endl;
  ofeff0 << "fraction ZERO/all : " << ineff_zero_S0.first << " +- " << ineff_zero_S0.second << endl;
  //****
  ofeff0 << "#presel_2mu_S1_passing_golden_events: " << N_presel_passing_2mu_golden_S1_events << endl;
  ofeff0 << "#2mu_S1_passing_golden_&_S0_trackable: " << N_presel_passing_2mu_golden_S1_trackable_S0_events << endl;
  ofeff0 << "#2mu_passing_golden_in_S0_&&_S1: " << N_passing_2mu_golden_events << endl;
  ofeff0 << "2mu efficiency S0 trackable / presel S1 passing golden = " << eff_presel_2mu_S0.first << " +- " << eff_presel_2mu_S0.second << endl;
  ofeff0 << "2mu efficiency S0 golden    / presel S1 pass.golden    = " << eff_golden_2mu_S0.first << " +- " << eff_golden_2mu_S0.second << endl;
  //
  ofeff0 << "#presel_3mu_S1_passing_golden_events: " << N_presel_passing_3mu_golden_S1_events << endl;
  ofeff0 << "#3mu_S1_passing_golden_&_S0_trackable: " << N_presel_passing_3mu_golden_S1_trackable_S0_events << endl;
  ofeff0 << "#3mu_passing_golden_in_S0_&&_S1: " << N_passing_3mu_golden_events << endl;
  ofeff0 << "3mu efficiency S0 trackable / presel S1 passing golden = " << eff_presel_3mu_S0.first << " +- " << eff_presel_3mu_S0.second << endl;
  ofeff0 << "3mu efficiency S0 golden    / presel S1 pass.golden    = " << eff_golden_3mu_S0.first << " +- " << eff_golden_3mu_S0.second << endl;
  //
  ofeff0 << "#presel_4mu_S1_passing_golden_events: " << N_presel_passing_4mu_golden_S1_events << endl;
  ofeff0 << "#4mu_S1_passing_golden_&_S0_trackable: " << N_presel_passing_4mu_golden_S1_trackable_S0_events << endl;
  ofeff0 << "#4mu_passing_golden_in_S0_&&_S1: " << N_passing_4mu_golden_events << endl;
  ofeff0 << "4mu efficiency S0 trackable / presel S1 passing golden = " << eff_presel_4mu_S0.first << " +- " << eff_presel_4mu_S0.second << endl;
  ofeff0 << "4mu efficiency S0 golden    / presel S1 pass.golden    = " << eff_golden_4mu_S0.first << " +- " << eff_golden_4mu_S0.second << endl;
  ofeff0.close();
  ////////////////////////////////////////////////////

  cout << endl;
  cout << "***** RATES w.r.t. all MERGED events ***************************************************" << endl;
  cout << "Fraction Passing 1mu             = " << double(N_passing_1mu_events) / double(N_input_events) << endl;
  cout << "Fraction Passing GOLDEN 1mu      = " << double(N_passing_1mu_golden_events) / double(N_input_events) << endl;
  cout << "Fraction Passing GOLDEN 2mu      = " << double(N_passing_2mu_golden_events) / double(N_input_events) << endl;
  cout << "Fraction Passing GOLDEN 3mu      = " << double(N_passing_3mu_golden_events) / double(N_input_events) << endl;
  cout << "Fraction Passing GOLDEN 4mu      = " << double(N_passing_4mu_golden_events) / double(N_input_events) << endl;
  cout << "Fraction Single Mu Int.          = " << double(N_single_clean_events) / double(N_input_events) << endl;
  cout << "Fraction Golden sel              = " << double(N_golden_events) / double(N_input_events) << endl;
  cout << "Fraction 2 Mu PU Int.            = " << double(N_pileup_2mu_events) / double(N_input_events) << endl;
  cout << "Fraction 3 Mu PU Int.            = " << double(N_pileup_3mu_events) / double(N_input_events) << endl;
  cout << "Fraction 4 Mu PU Int.            = " << double(N_pileup_4mu_events) / double(N_input_events) << endl;
  cout << "Fraction >=4 Mu PU Int.          = " << double(N_pileup_many_events) / double(N_input_events) << endl;
  cout << "Total Fraction PU 2+3+4 Mu Int.  = " << double(N_pileup234_skim_events) / double(N_input_events) << endl;
  cout << "Total Fraction PU Mu Int.        = " << double(N_pileup_skim_events) / double(N_input_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4   = " << double(N_loose234_events) / double(N_input_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4+  = " << double(N_loose_events) / double(N_input_events) << endl;
  cout << endl;
  cout << "===========================================================================================" << endl;
  cout << "Fraction of PRESELECTED events (with at least one stub in S0 and one stub in S1): "
       << double(N_presel_events) / double(N_input_events) << endl;
  cout << "===========================================================================================" << endl;
  cout << endl;
  cout << "***** RATES w.r.t. PRESELECTED events ***************************************************" << endl;
  cout << "Fraction Passing 1mu             = " << double(N_passing_1mu_events) / double(N_presel_events) << endl;
  cout << "Fraction Passing GOLDEN 1mu      = " << double(N_passing_1mu_golden_events) / double(N_presel_events) << endl;
  cout << "Fraction Passing GOLDEN 2mu      = " << double(N_passing_2mu_golden_events) / double(N_presel_events) << endl;
  cout << "Fraction Passing GOLDEN 3mu      = " << double(N_passing_3mu_golden_events) / double(N_presel_events) << endl;
  cout << "Fraction Passing GOLDEN 4mu      = " << double(N_passing_4mu_golden_events) / double(N_presel_events) << endl;
  cout << "Fraction Single Mu Int.          = " << double(N_single_clean_events) / double(N_presel_events) << endl;
  cout << "Fraction Golden sel              = " << double(N_golden_events) / double(N_presel_events) << endl;
  cout << "Fraction 2 Mu PU Int.            = " << double(N_pileup_2mu_events) / double(N_presel_events) << endl;
  cout << "Fraction 3 Mu PU Int.            = " << double(N_pileup_3mu_events) / double(N_presel_events) << endl;
  cout << "Fraction 4 Mu PU Int.            = " << double(N_pileup_4mu_events) / double(N_presel_events) << endl;
  cout << "Fraction >=4 Mu PU Int.          = " << double(N_pileup_many_events) / double(N_presel_events) << endl;
  cout << "Total Fraction PU 2+3+4 Mu Int.  = " << double(N_pileup234_skim_events) / double(N_presel_events) << endl;
  cout << "Total Fraction PU Mu Int.        = " << double(N_pileup_skim_events) / double(N_presel_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4   = " << double(N_loose234_events) / double(N_presel_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4+  = " << double(N_loose_events) / double(N_presel_events) << endl;
  cout << endl;
  cout << "===========================================================================================" << endl;
  cout << "Fraction of S0-Trackable events (PRESELECTED, with at least one stub in S1) w.r.t. all MERGED: "
       << double(N_trackable_0_onehit_1) / double(N_input_events) << endl;
  cout << "Fraction of S0-Trackable events (PRESELECTED, with at least one stub in S1) w.r.t. PRESELECTED: "
       << double(N_trackable_0_onehit_1) / double(N_presel_events) << endl;
  cout << "===========================================================================================" << endl;
  cout << endl;
  cout << "***** RATES w.r.t. S0-Trackable PRESELECTED events **************************************" << endl;
  cout << "Fraction Passing 1mu             = " << double(N_passing_1mu_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction Passing GOLDEN 1mu      = " << double(N_passing_1mu_golden_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction Passing GOLDEN 2mu      = " << double(N_passing_2mu_golden_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction Passing GOLDEN 3mu      = " << double(N_passing_3mu_golden_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction Passing GOLDEN 4mu      = " << double(N_passing_4mu_golden_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction Single Mu Int.          = " << double(N_single_clean_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction Golden sel              = " << double(N_golden_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction 2 Mu PU Int.            = " << double(N_pileup_2mu_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction 3 Mu PU Int.            = " << double(N_pileup_3mu_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction 4 Mu PU Int.            = " << double(N_pileup_4mu_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction >=4 Mu PU Int.          = " << double(N_pileup_many_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Total Fraction PU 2+3+4 Mu Int.  = " << double(N_pileup234_skim_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Total Fraction PU Mu Int.        = " << double(N_pileup_skim_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Total Fraction Single+PU 2+3+4   = " << double(N_loose234_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Total Fraction Single+PU 2+3+4+  = " << double(N_loose_events) / double(N_trackable_0_onehit_1) << endl;
  cout << endl;
  cout << "===========================================================================================" << endl;
  cout << "Fraction of TRACKABLE events w.r.t. all Merged events: "
       << double(N_trackable_events) / double(N_input_events) << endl;
  cout << "Fraction of TRACKABLE events w.r.t. PRESELECTED events: "
       << double(N_trackable_events) / double(N_presel_events) << endl;
  cout << "Fraction of TRACKABLE events w.r.t. PRESELECTED && S0-Trackable events: "
       << double(N_trackable_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "===========================================================================================" << endl;
  cout << endl;
  cout << "***** RATES w.r.t. TRACKABLE events ***************************************************" << endl;
  cout << "Fraction Passing 1mu             = " << double(N_passing_1mu_events) / double(N_trackable_events) << endl;
  cout << "Fraction Passing GOLDEN 1mu      = " << double(N_passing_1mu_golden_events) / double(N_trackable_events) << endl;
  cout << "Fraction Passing GOLDEN 2mu      = " << double(N_passing_2mu_golden_events) / double(N_trackable_events) << endl;
  cout << "Fraction Passing GOLDEN 3mu      = " << double(N_passing_3mu_golden_events) / double(N_trackable_events) << endl;
  cout << "Fraction Passing GOLDEN 4mu      = " << double(N_passing_4mu_golden_events) / double(N_trackable_events) << endl;
  cout << "Fraction Single Mu Int.          = " << double(N_single_clean_events) / double(N_trackable_events) << endl;
  cout << "Fraction Golden sel              = " << double(N_golden_events) / double(N_trackable_events) << endl;
  cout << "Fraction 2 Mu PU Int.            = " << double(N_pileup_2mu_events) / double(N_trackable_events) << endl;
  cout << "Fraction 3 Mu PU Int.            = " << double(N_pileup_3mu_events) / double(N_trackable_events) << endl;
  cout << "Fraction 4 Mu PU Int.            = " << double(N_pileup_4mu_events) / double(N_trackable_events) << endl;
  cout << "Fraction >=4 Mu PU Int.          = " << double(N_pileup_many_events) / double(N_trackable_events) << endl;
  cout << "Total Fraction PU 2+3+4 Mu Int.  = " << double(N_pileup234_skim_events) / double(N_trackable_events) << endl;
  cout << "Total Fraction PU Mu Int.        = " << double(N_pileup_skim_events) / double(N_trackable_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4   = " << double(N_loose234_events) / double(N_trackable_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4+  = " << double(N_loose_events) / double(N_trackable_events) << endl;

  cout << "===========================================================================================" << endl;
  cout << "Fraction of 1mu PASSING GOLDEN events w.r.t. all Merged events: "
       << double(N_passing_1mu_golden_events) / double(N_input_events) << endl;
  cout << "Fraction of 1mu PASSING GOLDEN events w.r.t. PRESELECTED events: "
       << double(N_passing_1mu_golden_events) / double(N_presel_events) << endl;
  cout << "Fraction of 1mu PASSING GOLDEN events w.r.t. PRESELECTED && S0-Trackable events: "
       << double(N_passing_1mu_golden_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction of 1mu PASSING GOLDEN events w.r.t. TRACKABLE events: "
       << double(N_passing_1mu_golden_events) / double(N_trackable_events) << endl;

  cout << "===========================================================================================" << endl;
  cout << endl;
  cout << "***** RATES w.r.t. 1mu PASSING GOLDEN events **********************************************" << endl;
  cout << "Fraction Passing 1mu             = " << double(N_passing_1mu_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction Passing GOLDEN 1mu      = " << double(N_passing_1mu_golden_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction Passing GOLDEN 2mu      = " << double(N_passing_2mu_golden_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction Passing GOLDEN 3mu      = " << double(N_passing_3mu_golden_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction Passing GOLDEN 4mu      = " << double(N_passing_4mu_golden_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction Single Mu Int.          = " << double(N_single_clean_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction Golden sel              = " << double(N_golden_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction 2 Mu PU Int.            = " << double(N_pileup_2mu_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction 3 Mu PU Int.            = " << double(N_pileup_3mu_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction 4 Mu PU Int.            = " << double(N_pileup_4mu_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Fraction >=4 Mu PU Int.          = " << double(N_pileup_many_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Total Fraction PU 2+3+4 Mu Int.  = " << double(N_pileup234_skim_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Total Fraction PU Mu Int.        = " << double(N_pileup_skim_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4   = " << double(N_loose234_events) / double(N_passing_1mu_golden_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4+  = " << double(N_loose_events) / double(N_passing_1mu_golden_events) << endl;

  cout << "===========================================================================================" << endl;
  cout << "Fraction of events with ALL 12 modules fired w.r.t. all Merged events: "
       << double(N_fired12_events) / double(N_input_events) << endl;
  cout << "Fraction of events with ALL 12 modules fired w.r.t. PRESELECTED events: "
       << double(N_fired12_events) / double(N_presel_events) << endl;
  cout << "Fraction of events with ALL 12 modules fired w.r.t. PRESELECTED && S0-Trackable events: "
       << double(N_fired12_events) / double(N_trackable_0_onehit_1) << endl;
  cout << "Fraction of events with ALL 12 modules fired w.r.t. TRACKABLE events: "
       << double(N_fired12_events) / double(N_trackable_events) << endl;

  cout << "===========================================================================================" << endl;
  cout << endl;
  cout << "***** RATES w.r.t. events with ALL 12 modules fired ***************************************" << endl;
  cout << "Fraction Passing 1mu             = " << double(N_passing_1mu_events) / double(N_fired12_events) << endl;
  cout << "Fraction Passing GOLDEN 1mu      = " << double(N_passing_1mu_golden_events) / double(N_fired12_events) << endl;
  cout << "Fraction Passing GOLDEN 2mu      = " << double(N_passing_2mu_golden_events) / double(N_fired12_events) << endl;
  cout << "Fraction Passing GOLDEN 3mu      = " << double(N_passing_3mu_golden_events) / double(N_fired12_events) << endl;
  cout << "Fraction Passing GOLDEN 4mu      = " << double(N_passing_4mu_golden_events) / double(N_fired12_events) << endl;
  cout << "Fraction Single Mu Int.          = " << double(N_single_clean_events) / double(N_fired12_events) << endl;
  cout << "Fraction Golden sel              = " << double(N_golden_events) / double(N_fired12_events) << endl;
  cout << "Fraction 2 Mu PU Int.            = " << double(N_pileup_2mu_events) / double(N_fired12_events) << endl;
  cout << "Fraction 3 Mu PU Int.            = " << double(N_pileup_3mu_events) / double(N_fired12_events) << endl;
  cout << "Fraction 4 Mu PU Int.            = " << double(N_pileup_4mu_events) / double(N_fired12_events) << endl;
  cout << "Fraction >=4 Mu PU Int.          = " << double(N_pileup_many_events) / double(N_fired12_events) << endl;
  cout << "Total Fraction PU 2+3+4 Mu Int.  = " << double(N_pileup234_skim_events) / double(N_fired12_events) << endl;
  cout << "Total Fraction PU Mu Int.        = " << double(N_pileup_skim_events) / double(N_fired12_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4   = " << double(N_loose234_events) / double(N_fired12_events) << endl;
  cout << "Total Fraction Single+PU 2+3+4+  = " << double(N_loose_events) / double(N_fired12_events) << endl;
  cout << "Fraction Single Mu Int.+FIRED12  = " << double(N_fired12_plus_single_clean_events) / double(N_fired12_events) << endl;

  // reset stat errors on the per module histos
  h_nstubsPerModule->Sumw2(false);
  h_nstubsPerModule_presel->Sumw2(false);
  h_nstubsPerModule_trackable_0_onehit_1->Sumw2(false);
  h_nstubsPerModule_trackable_1_onehit_0->Sumw2(false);
  h_nstubsPerModule_trackable->Sumw2(false);
  h_nstubsPerModule_fired12->Sumw2(false);
  h_nstubsPerModule_single_clean->Sumw2(false);
  h_nstubsPerModule_pileup234_skim->Sumw2(false);
  h_nstubsPerModule_fired12_plus_single_clean->Sumw2(false);

  hfile->Write();
    
  if (iskim > 0) {
    o_file = o_tree->GetCurrentFile();
    o_file->cd();
    o_tree->Write();
    o_file->Close();
    cout << "\n Number of output skimmed events          = " << o_nevt << endl;
    cout << "\n Number of stubs in output skimmed events = " << N_stubs_skim << endl;
  }

  return 0;
}
