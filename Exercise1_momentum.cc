#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TRandom.h"
#include "TMath.h"
#include "TF1.h"

double GetMomentum(double& grad_pre, double& grad_post) {
  double grad{0};
  if (grad_pre/grad_post <= 0) grad = abs(grad_pre-grad_post);
  else grad = abs(grad_pre+grad_post);
  double angle = TMath::ATan(grad);
  // p_y = 0.3 B * z * l , where B[T], z[e], l[m] and p[GeV/c]; p = p_y/sin(angle)
  // here B = 0.5T, l ~= 2m and z = 1e
  double p  =  0.3 / (2 * TMath::Sin(angle/2));
  return p; 
}

double GetGradientOfDCHits(std::vector<double> &branchX, std::vector<double> &branchZ, TRandom &rand) {
  //work in units of mm for x, m for z 
  double xerrors[5] {0.1, 0.1, 0.1, 0.1, 0.1};
  double zerrors[5] {0.0, 0.0, 0.0, 0.0, 0.0};
  double x[5] {0.0, 0.0, 0.0, 0.0, 0.0};
  double z[5] {0.0, 0.0, 0.0, 0.0, 0.0};

  
  for (Int_t j{0}; j < 5; j++) {
    x[j] = branchX[j] + rand.Gaus(0.0, 0.1);
    z[j] = branchZ[j]*0.5 + 2.75;
    std::cout << "Values for the point: " << x[j] << " and " << z[j] << std::endl;
    if (x[j] == 0.0 || z[j] == 0.0) std::cout << "WARNING: zero in coordinates" << std::endl;
  }
  /*
  ////////////////////////////////////////////
  // New DC hit selection procedure - method 1
  
  int old_z{-1};  //keep track of the previous z, controls indexing of kept values
  double old_x{0};   //kepp track of last saved x values
  int counter{0};    //index of entry in branch

  //iterate until all x-values are filled
  while (x[4] == 0.0) {
    
    //or until all hits for tree entry are iterated over
    if (branchZ.size() == counter-1) {
      //return a non-physical grad
      return 0.0;
    }
    
    //in the first DC layer, save the first hit for which x is non-zero
    else if (old_z == -1 && branchX[counter] != 0.0) {
      old_z += 1;
      x[old_z] = branchX[counter] + rand.Gaus(0.0, 0.1);
      z[old_z] = branchZ[counter] * 0.5 + 2.75;
      old_x = x[old_z];
    }
    
    //in subsequnt DC layers
    else if (branchZ[counter] == (double)old_z + 1.0) {
      //assign a hit if the current x is farther from zero than the last one saved
      if ( ((x[old_z] < 0.0) && (x[old_z] > branchX[counter])) || ((x[old_z] > 0.0) && (x[old_z] < branchX[counter])) ){
        old_z += 1;
        x[old_z] = branchX[counter] + rand.Gaus(0.0, 0.1);
        z[old_z] = branchZ[counter] * 0.5 + 2.75;
        old_x = x[old_z];
      }	
    }
    
    //move over to the next entry
    counter += 1;
  }
  */
  //fitting and extracting the gradient
  TGraphErrors* tge = new TGraphErrors(5, z, x, zerrors, xerrors);
  tge->Fit("pol1");
  return abs(tge->GetFunction("pol1")->GetParameter(1) / 1000); // convert the units

}

///////////////////////////////////////////////
// Track fitting and chi2 evaluation - Method 2

double GetGradientOfDCHitsCombo(std::vector<double> &branchX, std::vector<double> &branchZ, TRandom &rand) {
  //assign hits in each DC plane to vectors
  std::vector<double> x1, x2, x3, x4, x5;
  for (Int_t i{0}; i < branchZ.size(); i++) {
    if (branchZ[i] == 0)
      x1.push_back(branchX[i] + rand.Gaus(0.0, 0.1));
    else if (branchZ[i] == 1)
      x2.push_back(branchX[i] + rand.Gaus(0.0, 0.1));
    else if (branchZ[i] == 2)
      x3.push_back(branchX[i] + rand.Gaus(0.0, 0.1));
    else if (branchZ[i] == 3)
      x4.push_back(branchX[i] + rand.Gaus(0.0, 0.1));
    else if (branchZ[i] == 4)
      x5.push_back(branchX[i] + rand.Gaus(0.0, 0.1));
  }
  
  //vars to keep track of best chi2 and correposnding gradient value
  double c2_min = 100, best_m = 0;
  //constant values for fits
  double xerrors[5] {0.1, 0.1, 0.1, 0.1, 0.1};
  double zerrors[5] {0.0, 0.0, 0.0, 0.0, 0.0};
  double z[5] {2.75, 3.25, 3.75, 4.25, 4.75};

  //iterate over all possible combinations
  for (int i{0}; i < x1.size(); i++){
    for (int j{0}; j < x2.size(); j++){
      for (int k{0}; k < x3.size(); k++){
        for (int l{0}; l < x4.size(); l++){
          for (int m{0}; m < x5.size(); m++){
            double x[5] {x1[i], x2[j], x3[k], x4[l], x5[m]};  //set x values
            TGraphErrors* tge = new TGraphErrors(5, z, x, zerrors, xerrors);
            tge->Fit("pol1");
	    //save values if this is the best chi2 ever found
            if (tge->GetFunction("pol1")->GetChisquare() < c2_min) {
	      c2_min = tge->GetFunction("pol1")->GetChisquare();
	      best_m = tge->GetFunction("pol1")->GetParameter(1) / 1000; // convert the units
	    }
	  }
	}
      }
    }
  }
  
  // if all fits are poor, retrun 0
  if (c2_min < 10) return best_m;
  else return 0.0;

}

bool MatchTrackToHodCol(TF1 &fit, int& colNo){
  // ALLOW FOR DIFFERENT SIGNS OF M and C
  // try to return a non-zero number (distance to tile if not exact hit)
  // to give a measure of closeness of candidate
  double z_hod{5.0};
  double grad =  fit.GetParameter(1);
  double intercept = fit.GetParameter(0);
  // exactly
  if ( ((int)(z_hod * grad + intercept))/100 == colNo)
    return true;
  // +1 sigma for grad and intercept
  else if ( ((int)(z_hod * (grad+fit.GetParError(1)) + intercept+fit.GetParError(0)))/100 == colNo)
    return true;
  // -1 sigma for grad and intercept
  else if ( ((int)(z_hod * (grad-fit.GetParError(1)) + intercept-fit.GetParError(0)))/100 == colNo)
    return true;
  // no match
  else return false;
}

///////////////////////////////////
// Main function

void Exercise1_momentum() {
  
  TFile* b5 = new TFile("./B5_THIN_PbPlate_after_1000mus_100GeV_05T_0deg.root", "read");
  TTree* t = (TTree*) b5->Get("B5");
  std::vector<double> *HitXpre = new std::vector<double>;
  std::vector<double> *HitZpre = new std::vector<double>;
  std::vector<double> *HitX = new std::vector<double>;
  std::vector<double> *HitZ = new std::vector<double>;
  t->SetBranchAddress("Dc1HitsVector_x", &HitXpre);  
  t->SetBranchAddress("Dc1HitsVector_z", &HitZpre);
  t->SetBranchAddress("Dc2HitsVector_x", &HitX);  
  t->SetBranchAddress("Dc2HitsVector_z", &HitZ);
 
  TRandom* rand = new TRandom();

  TH1F* histo = new TH1F("THIN PbPlate after B=0.5T", "Momentum distribution B=0.5T with THIN PbPLate", 50, 0.5, 150.5); 
  histo->GetXaxis()->SetTitle("Momentum [GeV]");
  histo->GetYaxis()->SetTitle("Count per 3.0 GeV bin"); 

  int unsucessfull_hit_associations{0};
  int nonzero_DC1_gradient{0};
  std::vector<double> DC1_grads;  
  Long64_t Entries = t->GetEntries();
  for (Long64_t i{0}; i < Entries; i++) {
    t->GetEntry(i);
    double grad_pre = GetGradientOfDCHits(*HitXpre, *HitZpre, *rand);
    double grad_post = GetGradientOfDCHits(*HitX, *HitZ, *rand);    
    if (grad_post == 0.0) {
      unsucessfull_hit_associations += 1;
      continue;
    }
    else if (grad_pre != 0.0) {
      nonzero_DC1_gradient+=1;
      DC1_grads.push_back(grad_pre);
    }
    grad_pre = 0.0;
    double p = GetMomentum(grad_pre, grad_post);
    histo->Fill(p);
  }

  TFile* out = new TFile("Momentum_distribution_THIN_PbPlate_after_05T_0deg_2grads.root", "update");
  histo->Write();
  out->Close();
  delete histo;
  std::cout << "WARNING: " << unsucessfull_hit_associations << " out of " << Entries << " were unsuccessfull" << std::endl;
  std::cout << "with " << nonzero_DC1_gradient << " giving a non-zero gradient of tracks in the first DC of values: " << std::endl;
  /*  for (int i{0}; i < DC1_grads.size(); i++){
    std::cout << DC1_grads[i] << " /1000" << std::endl;
  }
  */
}
