#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"


#include "releff.h"

using namespace std;

int main(int argc, char *argv[])
{

    TString hname;

    if (argc != 2)
    {
        cout << "ERROR! Usage: " << argv[0] << " <hfile> " << endl;
        return 10;
    }
    else
    {
        hname = argv[1];
    }

    TFile *f_all = new TFile(hname + ".root");

    //  TH1D *h_nstubsPerModule = (TH1D*)f_all->Get("h_nstubsPerModule");
    TH1D *h_nstubsPerModule_S0 = (TH1D *)f_all->Get("h_nstubsPerModule_S0");
    TH1D *h_nstubsPerModule_S1 = (TH1D *)f_all->Get("h_nstubsPerModule_S1");

    TH1D *h_nTotStubs = (TH1D *)f_all->Get("h_nTotStubs");
    //  TH2D *h2_nStubs = (TH2D*)f_all->Get("h2_nStubs");
    //  TH1D *h_nStubs_0 = (TH1D*)f_all->Get("h_nStubs_0");
    //  TH1D *h_nStubs_1 = (TH1D*)f_all->Get("h_nStubs_1");

    TH1D *h_nTotStubs_presel = (TH1D *)f_all->Get("h_nTotStubs_presel");
    //  TH2D *h2_nStubs_presel = (TH2D*)f_all->Get("h2_nStubs_presel");
    //  TH1D *h_nStubs_0_presel = (TH1D*)f_all->Get("h_nStubs_0_presel");
    //  TH1D *h_nStubs_1_presel = (TH1D*)f_all->Get("h_nStubs_1_presel");

    TH1D *h_nTotStubs_S0_trackable_S1_onehit = (TH1D *)f_all->Get("h_nTotStubs_S0_trackable_S1_onehit");
    //  TH2D *h2_nStubs_S0_trackable_S1_onehit = (TH2D*)f_all->Get("h2_nStubs_S0_trackable_S1_onehit");
    //  TH1D *h_nStubs_0_S0_trackable_S1_onehit = (TH1D*)f_all->Get("h_nStubs_0_S0_trackable_S1_onehit");
    //  TH1D *h_nStubs_1_S0_trackable_S1_onehit = (TH1D*)f_all->Get("h_nStubs_1_S0_trackable_S1_onehit");

    TH1D *h_nTotStubs_trackable = (TH1D *)f_all->Get("h_nTotStubs_trackable");
    //  TH2D *h2_nStubs_trackable = (TH2D*)f_all->Get("h2_nStubs_trackable");
    //  TH1D *h_nStubs_0_trackable = (TH1D*)f_all->Get("h_nStubs_0_trackable");
    //  TH1D *h_nStubs_1_trackable = (TH1D*)f_all->Get("h_nStubs_1_trackable");

    //  TH1D *h_nTotStubs_zero_0 = (TH1D*)f_all->Get("h_nTotStubs_zero_0");
    //  TH1D *h_nTotStubs_zero_1 = (TH1D*)f_all->Get("h_nTotStubs_zero_1");
    TH1D *h_nTotStubs_noise = (TH1D *)f_all->Get("h_nTotStubs_noise");
    TH1D *h_nTotStubs_missing_0 = (TH1D *)f_all->Get("h_nTotStubs_missing_0");
    TH1D *h_nTotStubs_missing_1 = (TH1D *)f_all->Get("h_nTotStubs_missing_1");
    TH1D *h_nTotStubs_not_trackable = (TH1D *)f_all->Get("h_nTotStubs_not_trackable");
    //  TH1D *h_nTotStubs_presel_not_trackable = (TH1D*)f_all->Get("h_nTotStubs_presel_not_trackable");

    TH1D *h_nTotStubs_passing_1mu = (TH1D *)f_all->Get("h_nTotStubs_passing_1mu");

    TH1D *h_nTotStubs_golden = (TH1D *)f_all->Get("h_nTotStubs_golden");
    //  TH1D* h_nTotStubs_umsel = (TH1D*)f_all->Get("h_nTotStubs_umsel");
    TH1D *h_nTotStubs_single_cand = (TH1D *)f_all->Get("h_nTotStubs_single_cand");
    TH1D *h_nTotStubs_single_clean = (TH1D *)f_all->Get("h_nTotStubs_single_clean");
    //  TH2D* h2_nStubs_single_clean = (TH2D*)f_all->Get("h2_nStubs_single_clean");
    //  TH1D* h_nStubs_0_single_clean = (TH1D*)f_all->Get("h_nStubs_0_single_clean");
    //  TH1D* h_nStubs_1_single_clean = (TH1D*)f_all->Get("h_nStubs_1_single_clean");

    TH1D *h_nTotStubs_pileup_any = (TH1D *)f_all->Get("h_nTotStubs_pileup_any");
    TH1D *h_nTotStubs_pileup_2mu = (TH1D *)f_all->Get("h_nTotStubs_pileup_2mu");
    TH1D *h_nTotStubs_pileup_3mu = (TH1D *)f_all->Get("h_nTotStubs_pileup_3mu");
    TH1D *h_nTotStubs_pileup_4mu = (TH1D *)f_all->Get("h_nTotStubs_pileup_4mu");
    TH1D *h_nTotStubs_pileup_many = (TH1D *)f_all->Get("h_nTotStubs_pileup_many");
    TH1D *h_nTotStubs_pileup_skim = (TH1D *)f_all->Get("h_nTotStubs_pileup_skim");
    //  TH2D* h2_nStubs_pileup_skim = (TH2D*)f_all->Get("h2_nStubs_pileup_skim");
    //  TH1D* h_nStubs_0_pileup_skim = (TH1D*)f_all->Get("h_nStubs_0_pileup_skim");
    //  TH1D* h_nStubs_1_pileup_skim = (TH1D*)f_all->Get("h_nStubs_1_pileup_skim");
    TH1D *h_nTotStubs_pileup234_skim = (TH1D *)f_all->Get("h_nTotStubs_pileup234_skim");
    //  TH2D *h2_nStubs_pileup234_skim = (TH2D*)f_all->Get("h2_nStubs_pileup234_skim");
    //  TH1D *h_nStubs_0_pileup234_skim = (TH1D*)f_all->Get("h_nStubs_0_pileup234_skim");
    //  TH1D *h_nStubs_1_pileup234_skim = (TH1D*)f_all->Get("h_nStubs_1_pileup234_skim");
    TH1D *h_nTotStubs_fired12 = (TH1D *)f_all->Get("h_nTotStubs_fired12");
    //  TH2D *h2_nStubs_fired12 = (TH2D*)f_all->Get("h_Stubs_fired12");
    //  TH1D *h_nStubs_0_fired12 = (TH1D*)f_all->Get("h_Stubs_0_fired12");
    //  TH1D *h_nStubs_1_fired12 = (TH1D*)f_all->Get("h_Stubs_1_fired12");
    TH1D *h_nTotStubs_fired12_plus_single_clean = (TH1D *)f_all->Get("h_nTotStubs_fired12_plus_single_clean");
    //  TH2D *h2_nStubs_fired12_plus_single_clean = (TH2D*)f_all->Get("h2_nStubs_fired12_plus_single_clean");
    //  TH1D *h_nStubs_0_fired12_plus_single_clean = (TH1D*)f_all->Get("h_nStubs_0_fired12_plus_single_clean");
    //  TH1D *h_nStubs_1_fired12_plus_single_clean = (TH1D*)f_all->Get("h_nStubs_1_fired12_plus_single_clean");
    TH1D *h_passing_1mu_golden_vs_time = (TH1D *)f_all->Get("h_passing_1mu_golden_vs_time");
    TH1D *h_passing_2mu_golden_vs_time = (TH1D *)f_all->Get("h_passing_2mu_golden_vs_time");
    TH1D *h_passing_3mu_golden_vs_time = (TH1D *)f_all->Get("h_passing_3mu_golden_vs_time");
    TH1D *h_passing_4mu_golden_vs_time = (TH1D *)f_all->Get("h_passing_4mu_golden_vs_time");

    Long64_t N_input_events = h_nTotStubs->GetEntries();
    Long64_t N_input_stubs_0 = h_nstubsPerModule_S0->GetEntries();
    Long64_t N_input_stubs_1 = h_nstubsPerModule_S1->GetEntries();

    //  Long64_t N_stubs_skim = 0;
    Long64_t N_noise_events = h_nTotStubs_noise->GetEntries();
    Long64_t N_missing_0_events = h_nTotStubs_missing_0->GetEntries();
    Long64_t N_missing_1_events = h_nTotStubs_missing_1->GetEntries();
    Long64_t N_presel_events = h_nTotStubs_presel->GetEntries();
    Long64_t N_trackable_0_onehit_1 = h_nTotStubs_S0_trackable_S1_onehit->GetEntries();
    Long64_t N_not_trackable_events = h_nTotStubs_not_trackable->GetEntries();
    Long64_t N_trackable_events = h_nTotStubs_trackable->GetEntries();
    Long64_t N_fired12_events = h_nTotStubs_fired12->GetEntries();
    Long64_t N_passing_1mu_events = h_nTotStubs_passing_1mu->GetEntries();
    Long64_t N_passing_1mu_golden_events = h_passing_1mu_golden_vs_time->GetEntries();
    Long64_t N_passing_2mu_golden_events = h_passing_2mu_golden_vs_time->GetEntries();
    Long64_t N_passing_3mu_golden_events = h_passing_3mu_golden_vs_time->GetEntries();
    Long64_t N_passing_4mu_golden_events = h_passing_4mu_golden_vs_time->GetEntries();
    Long64_t N_single_cand_events = h_nTotStubs_single_cand->GetEntries();
    Long64_t N_single_clean_events = h_nTotStubs_single_clean->GetEntries();
    Long64_t N_golden_events = h_nTotStubs_golden->GetEntries();
    Long64_t N_pileup_any_events = h_nTotStubs_pileup_any->GetEntries();
    Long64_t N_pileup_2mu_events = h_nTotStubs_pileup_2mu->GetEntries();
    Long64_t N_pileup_3mu_events = h_nTotStubs_pileup_3mu->GetEntries();
    Long64_t N_pileup_4mu_events = h_nTotStubs_pileup_4mu->GetEntries();
    Long64_t N_pileup_many_events = h_nTotStubs_pileup_many->GetEntries();
    Long64_t N_pileup_skim_events = h_nTotStubs_pileup_skim->GetEntries();
    Long64_t N_pileup234_skim_events = h_nTotStubs_pileup234_skim->GetEntries();
    Long64_t N_loose_events = N_single_clean_events + N_pileup_skim_events;
    Long64_t N_loose234_events = N_single_clean_events + N_pileup234_skim_events;
    Long64_t N_fired12_plus_single_clean_events = h_nTotStubs_fired12_plus_single_clean->GetEntries();

    ifstream sumold1("relative_efficiency_S1_goodS0_sumold.txt");
    ifstream i1("relative_efficiency_S1_goodS0.txt");
    ofstream sum1("relative_efficiency_S1_goodS0_sum.txt");
    add_efficiency('1', '0', sumold1, i1, sum1);
    //  sum1.close();
    ifstream sumold0("relative_efficiency_S0_goodS1_sumold.txt");
    ifstream i0("relative_efficiency_S0_goodS1.txt");
    ofstream sum0("relative_efficiency_S0_goodS1_sum.txt");
    add_efficiency('0', '1', sumold0, i0, sum0);
    //  sum0.close();

    cout << "\n EVENTS:   " << endl;
    cout << "Number all input events = " << endl;
    cout << "\t" << N_input_events << endl;
    //  cout<<"Number of events without any hit in both S_0 and S_1 = "<<endl;
    //  cout<<"\t" << N_zero_events <<endl;
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

    cout << "Number of trackable events = " << endl; // both stations trackable
    cout << "\t " << N_trackable_events << endl;

    cout << "Number of events with ALL 12 fired modules = " << endl;
    cout << "\t" << N_fired12_events << endl;

    cout << "Number of clean 1mu passing events = " << endl;
    cout << "\t " << N_passing_1mu_events << endl;
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
    cout << "N input stubs tot = " << N_input_stubs_0 + N_input_stubs_1 << endl;

    cout << endl;
    cout << "Fraction missing S_0             = " << double(N_missing_0_events) / double(N_input_events) << endl;
    cout << "Fraction missing S_1             = " << double(N_missing_1_events) / double(N_input_events) << endl;
    cout << "Fraction not fully trackable     = " << double(N_not_trackable_events) / double(N_input_events) << endl;

    cout << endl
         << "N PRESELECTED Events     = " << N_presel_events << endl;

    cout << endl
         << "N PRESELECTED S_0 Trackable Events     = " << N_trackable_0_onehit_1 << endl;

    cout << endl
         << "N TRACKABLE Events       = " << N_trackable_events << endl;

    cout << "Fraction PRESELECTED / all       = " << double(N_presel_events) / double(N_input_events) << endl;
    cout << "Fraction PRESEL. S_0-Trackable / PRESEL = " << double(N_trackable_0_onehit_1) / double(N_presel_events) << endl;
    cout << "Fraction TRACKABLE / PRESEL      = " << double(N_trackable_events) / double(N_presel_events) << endl;

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

    //******

    TCanvas c;
    c.SetLogy();
    h_nTotStubs_trackable->SetTitle("");
    h_nTotStubs_trackable->GetXaxis()->SetTitle("Number of hits");
    h_nTotStubs_trackable->SetStats(0);
    h_nTotStubs_trackable->SetLineColor(1);
    h_nTotStubs_trackable->Draw();
    h_nTotStubs_single_clean->SetStats(0);
    h_nTotStubs_single_clean->SetLineColor(2);
    h_nTotStubs_single_clean->Draw("same");
    h_nTotStubs_pileup_2mu->SetStats(0);
    h_nTotStubs_pileup_2mu->SetLineColor(3);
    h_nTotStubs_pileup_2mu->Draw("same");
    h_nTotStubs_pileup_3mu->SetStats(0);
    h_nTotStubs_pileup_3mu->SetLineColor(4);
    h_nTotStubs_pileup_3mu->Draw("same");
    h_nTotStubs_pileup_4mu->SetStats(0);
    h_nTotStubs_pileup_4mu->SetLineColor(7);
    h_nTotStubs_pileup_4mu->Draw("same");

    TLegend *legend = new TLegend(0.50, 0.63, 0.88, 0.88);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_nTotStubs_trackable ,"Trackable events", "l");
    legend->AddEntry( h_nTotStubs_single_clean,"Single mu int", "l");
    legend->AddEntry(h_nTotStubs_pileup_2mu ,"2 mu pileup int", "l");
    legend->AddEntry(h_nTotStubs_pileup_3mu ,"3 mu pileup int", "l");
    legend->AddEntry(h_nTotStubs_pileup_4mu ,"4 mu pileup int", "l");

    // cAngBS->BuildLegend();
    legend->Draw();
    // c.BuildLegend();
    c.Modified();
    h_nTotStubs->SetTitle("");

    c.SaveAs("plot_nTotStubs.root");
    c.SaveAs("plot_nTotStubs.C");
    c.SaveAs("plot_nTotStubs.pdf");
    c.SaveAs("plot_nTotStubs.png");

    return 0;
}
