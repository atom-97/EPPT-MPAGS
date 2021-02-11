#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TRandom.h"
#include "TMath.h"
#include "TF1.h"
#include <algorithm>
#include <numeric>

void Exercise2_calorimeters() {
  
  TFile* b5 = new TFile("./B5_1000ps_100GeV_05T_0deg.root", "read");
  TTree* t = (TTree*) b5->Get("B5");
  std::vector<double> *EMCalVec = new std::vector<double>;
  std::vector<double> *HCalVec = new std::vector<double>;
  double EMCal;
  double HCal;
  t->SetBranchAddress("ECEnergyVector", &EMCalVec);  
  t->SetBranchAddress("HCEnergyVector", &HCalVec);
  t->SetBranchAddress("ECEnergy", &EMCal);  
  t->SetBranchAddress("HCEnergy", &HCal);

  TH1F* histo1 = new TH1F("EMtot p+ at 100GeV, B=0.5T", "Total energy deposit in EM Calorimeter by 100 GeV p+", 141, -0.25, 70.25); 
  histo1->GetXaxis()->SetTitle("Energy [GeV]");
  histo1->GetYaxis()->SetTitle("Count per 0.5 GeV bin"); 

  TH1F* histo2 = new TH1F("Htot p+ at 100GeV, B=0.5T", "Total energy deposit in Hadronic Calorimeter by 100 GeV p+", 201, -0.015, 6.015); 
  histo2->GetXaxis()->SetTitle("Energy [GeV]");
  histo2->GetYaxis()->SetTitle("Count per 0.03 GeV bin");

  TH1F* histo3 = new TH1F("EMavg p+ at 100GeV, B=0.5T", "Average energy deposit in EM Calorimeter by 100 GeV p+", 49, -0.25, 24.25); 
  histo3->GetXaxis()->SetTitle("Energy [GeV]");
  histo3->GetYaxis()->SetTitle("Count per 0.5 GeV bin"); 

  TH1F* histo4 = new TH1F("Havg p+ at 100GeV, B=0.5T", "Average energy deposit in Hadronic Calorimeter by 100 GeV p+", 81, -0.015, 2.415); 
  histo4->GetXaxis()->SetTitle("Energy [GeV]");
  histo4->GetYaxis()->SetTitle("Count per 0.03 GeV bin");

  TH1F* histo5 = new TH1F("EMn p+ at 100GeV, B=0.5T", "Number of hit cells in EM Calorimeter by 100 GeV p+", 81, -0.5, 80.5); 
  histo5->GetXaxis()->SetTitle("Number of hit cells per incident particle [1]");
  histo5->GetYaxis()->SetTitle("Count");

  TH1F* histo6 = new TH1F("Hn p+ at 100GeV, B=0.5T", "Number of hit cells in Hadronic Calorimeter by 100 GeV p+", 21, -0.5, 20.5); 
  histo6->GetXaxis()->SetTitle("Number of hit cells per incident particle [1]");
  histo6->GetYaxis()->SetTitle("Count");

  TH2F* histo7 = new TH2F("EMnvsavg p+ at 100GeV, B=0.5T", "Number of hit cells vs average energy deposit in EM Calorimeter by 100 GeV p+", 49, -0.25, 24.25, 81, -0.5, 80.5);
  histo7->GetYaxis()->SetTitle("Number of hit cells per incident particle [1]");
  histo7->GetXaxis()->SetTitle("Energy [GeV per 0.5 GeV bin]");

  TH2F* histo8 = new TH2F("Hnvsavg p+ at 100GeV, B=0.5T", "Number of hit cells vs average energy deposit in Hadronic Calorimeter by 100 GeV p+", 81, -0.015, 2.415, 21, -0.5, 20.5); 
  histo8->GetYaxis()->SetTitle("Number of hit cells per incident particle [1]");
  histo8->GetXaxis()->SetTitle("Energy [GeV per 0.03 GeV bin]");
  
  int allEMCellsHit{0};
  int allHCellsHit{0};
  Long64_t Entries = t->GetEntries();
  for (Long64_t i{0}; i < Entries; i++) {
    t->GetEntry(i);

    //all values stored in GeV for mu+
    //EM in GeV and H in MeV for p+
    //all in GeV for p+
    //hence the conversion factor of 1/1000
    histo1->Fill(EMCal/1000);
    histo2->Fill(HCal/1000);
    auto nEM = EMCalVec->size() - std::count(EMCalVec->begin(), EMCalVec->end(), 0.0);
    if (nEM == EMCalVec->size()) allEMCellsHit += 1; 
    auto nHad = HCalVec->size() - std::count(HCalVec->begin(), HCalVec->end(), 0.0);
    if (nHad == HCalVec->size()) allHCellsHit += 1;
    histo3->Fill(EMCal/nEM/1000);
    histo4->Fill(HCal/nHad/1000);
    histo5->Fill(nEM);
    histo6->Fill(nHad);
    histo7->Fill(EMCal/nEM/1000,nEM);
    histo8->Fill(HCal/nHad/1000,nHad);
  }
  TFile* out = new TFile("./Cal_Energy_Distros_1000p+.root", "update");
  histo1->Write();
  histo2->Write();
  histo3->Write();
  histo4->Write();
  histo5->Write();
  histo6->Write();
  histo7->Write();
  histo8->Write();
  out->Close();
  std::cout << "Note:" << std::endl << " in " << allEMCellsHit << " out of " << Entries << " cases all cellss in the EM cal had energy deposits " << std::endl << " and in " << allHCellsHit << " out of " << Entries << " cases the same was true for Had cal" << std::endl; 
}
