#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLine.h"
#include "TF1.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TFeldmanCousins.h"
#include <TApplication.h>

using namespace std;

const int NMAX = 100000;

int main( int argc, char* argv[]){
  //TApplication * myapp = new TApplication("myapp",0,0);
  char *resdir = argv[1];
  char *tier_name = argv[2];
  int chn = atoi(argv[3]);
  bool dplms = true;
  if (argc>4) dplms = atoi(argv[4]);
  
  cout << "Using file " << tier_name << endl;
  cout << "Analizing channel " << chn << endl;
  if (dplms) cout << "DPLMS analysis " << endl;
  else cout << "Standard analysis " << endl;
  
  TFile *tier_file = TFile::Open(tier_name);
  TTree *tier = (TTree*)tier_file->Get("tier");
  
  std::vector<double>* energy = NULL;    
  tier->SetBranchAddress("energy", &energy);    
  std::vector<double>* AoE = NULL;
  tier->SetBranchAddress("AoE", &AoE);    
  std::vector<double>* AoE_dplms = NULL;
  tier->SetBranchAddress("AoE_dplms", &AoE_dplms);    
  ULong64_t tempo;
  tier->SetBranchAddress ("tempo",  &tempo);
  
  double mine = 1;
  double maxe = 3000;
  int bin_in = 10000;
  std::vector<double>  energies;
  std::vector<double>  AoEs;
  std::vector<double>  times;
  
  const int nEntries = tier->GetEntries();
  cout << "Total events = " << nEntries << endl;
  clock_t begin, now, time;
  begin = clock();
  //for ( int iEntry = 0; iEntry < 5000000;iEntry++){
  for ( int iEntry = 0; iEntry < nEntries; iEntry++ ) {
    if ((iEntry+1) % 1000000 == 0) cout << "event n." << iEntry+1 << endl;
    tier->GetEntry(iEntry);
    if ( energy->at(chn) > 900 ){
      energies.push_back(energy->at(chn));
      times.push_back(tempo);
      if ( dplms ) AoEs.push_back(AoE_dplms->at(chn));
      else AoEs.push_back(AoE->at(chn));
    }
  }
  int events = energies.size();
  cout << "Selected events: " << events << endl;
  now = clock();
  time = now - begin;
  cout << "The time to read waveforms: " << time/1.e6 << " s" << endl;
  
  TH1D *hene = new TH1D("hene","; Energy [keV]; Counts", bin_in, mine, maxe);
  for ( int i = 0; i < events;i++)
    hene->Fill(energies.at(i));

  //---------------- Print energy spectrum -----------------------  
  
  TCanvas *c = new TCanvas("c","energy");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0); 
  gPad->SetLogy();
  hene->Draw();
  c->Update();
  c->Print(Form("%s/chn%d_spectrum.pdf",resdir,chn));    
  
  
  //----- fit of the energy peaks with a double gaussian --------
  double min = 1570;
  double max = 1640;
  TF1 *f_2gaus = new TF1("f_2gaus","gaus(0) + (1. - 0.5*TMath::Erfc((x - [1])/[2]/TMath::Sqrt(2.0))) + gaus(3)+pol0(6)",min ,max);  
  
  // initialization of parameters of the fit
  hene->GetXaxis()->SetRangeUser(min,(max+min)/2);
  double A_d = hene->GetMaximum();
  double mu_d = hene->GetBinCenter( hene->GetMaximumBin());
  hene->GetXaxis()->SetRangeUser((max+min)/2,max);
  double A_f = hene->GetMaximum();
  double mu_f = hene->GetBinCenter( hene->GetMaximumBin());
    
  f_2gaus->SetParameter(0, A_d);
  f_2gaus->SetParameter(1, mu_d);
  f_2gaus->SetParameter(2, mu_d * 0.001);
  f_2gaus->SetParLimits(2, 0, 100);
  f_2gaus->SetParameter(3, A_f);
  f_2gaus->SetParameter(4, mu_f);
  f_2gaus->SetParameter(5, mu_f * 0.001);
  f_2gaus->SetParLimits(5, 0, 100);
  
  hene->GetXaxis()->SetRangeUser(1570,1585);
  double max_bkg = hene->GetMaximum();
  double min_bkg = hene->GetMinimum();
  double bkg_ini = (max_bkg+min_bkg)/2.;
  f_2gaus->SetParameter(6, bkg_ini);
  
  cout << "Initialization of DEP&FEP background  = " << bkg_ini << endl;
  hene->GetXaxis()->SetRangeUser(min,max);  
  gStyle->SetOptFit(0);
  hene->Fit("f_2gaus", "L", "", min, max);
  c->Update();
  c->Print(Form("%s/chn%d_DEP_FEP.pdf",resdir,chn));
  
  double A_dep = f_2gaus->GetParameter(0);
  double mu_dep = f_2gaus->GetParameter(1);
  double sigma_dep = f_2gaus->GetParameter(2);
  double A_fep = f_2gaus->GetParameter(3);
  double mu_fep = f_2gaus->GetParameter(4);
  double sigma_fep = f_2gaus->GetParameter(5);
  double bkg = f_2gaus->GetParameter(6);
  
  double n_sigma = 4.5;//number of sigma to calculate area  

  //---- DEP interval ----
  double min_dep = mu_dep - n_sigma * sigma_dep;
  double max_dep = mu_dep + n_sigma * sigma_dep;
  double min_dep2 = mu_dep - 2*n_sigma * sigma_dep;
  double max_dep2 = mu_dep + 2*n_sigma * sigma_dep;
  double max_dep3 = mu_dep + 3*n_sigma * sigma_dep;

  ////////////////////////////////////////////////////
  //----------------- AoE_analysis ------------------
  ////////////////////////////////////////////////////
  
  int countdep = 0;
  TH1D *hAoE0 = new TH1D("hAoE0","; A/E; Counts",1000,0,0.2);  
  for ( int i = 0; i < events;i++){
    if ( (energies.at(i) >= min_dep) && (energies.at(i) <= max_dep))
      countdep++;
      hAoE0->Fill(AoEs.at(i));
  }
  cout << "DEP events " << countdep << endl;
  TCanvas *c1 = new TCanvas("c1","AoE");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);
  hAoE0->Draw();
  c1->Update();
  c1->Print(Form("%s/chn%d_AoE_all.pdf",resdir,chn));
  
  double mean_aoe = hAoE0->GetBinCenter(hAoE0->GetMaximumBin());
  cout << "mean_aoe = " << mean_aoe << endl;
  
  //--------------- draw the A/E of the DEP --------------
  int nbin = 100;
  int nbin_cut = 1000;
  TH1D *hAoE = new TH1D("hAoE",Form("chn%d;A/E;Counts",chn), nbin, 0.7, 1.1);
  TH1D *hAoE1 = new TH1D("hAoE1",Form(";A/E;Counts",chn), nbin, 0.7, 1.1);
  TH1D *hAoE2 = new TH1D("hAoE2",Form(";A/E;Counts",chn), nbin, 0.7, 1.1);
  TH1D *hAoE3 = new TH1D("hAoE3", Form(";A/E;Counts",chn), nbin, 0.7, 1.1);
  TH1D *h_AoE = new TH1D("h_AoE",Form("chn%d;A/E;Counts",chn), nbin_cut, 0.7, 1.1);
  TH1D *h_AoE1 = new TH1D("h_AoE1",Form(";A/E;Counts",chn), nbin_cut, 0.7, 1.1);
  TH1D *h_AoE2 = new TH1D("h_AoE2",Form(";A/E;Counts",chn), nbin_cut, 0.7, 1.1);
  TH1D *h_AoE3 = new TH1D("h_AoE3", Form(";A/E;Counts",chn), nbin_cut, 0.7, 1.1);
  for ( int i = 0; i < events;i++){
    double AoEnorm = AoEs.at(i)/mean_aoe;
    if ( (energies.at(i) >= min_dep) && (energies.at(i) <= max_dep) ) {
      hAoE->Fill(AoEnorm);
      h_AoE->Fill(AoEnorm);
    }
    else if ( (energies.at(i) >= min_dep2) && (energies.at(i) <= min_dep) ) {
      hAoE1->Fill(AoEnorm);
      h_AoE1->Fill(AoEnorm);
    }
    else if ( (energies.at(i) >= max_dep) && (energies.at(i) <= max_dep2) ) {
      hAoE2->Fill(AoEnorm);
      h_AoE2->Fill(AoEnorm);
    }
    else if ( (energies.at(i) >= max_dep2) && (energies.at(i) <= max_dep3) ) {
      hAoE3->Fill(AoEnorm);
      h_AoE3->Fill(AoEnorm);
    }
  }
  
  hAoE2->Add(hAoE2);
  hAoE1->Add(hAoE2);
  hAoE1->Add(hAoE3,-1);
  hAoE->Add(hAoE1,-1);
  hAoE->Scale(1./hAoE->GetBinContent(hAoE->GetMaximumBin()));
  hAoE->SetLineWidth(2);
  hAoE->Draw();
  c1->Update();
  
  //---------------- fit the A/E of the DEP ------------
  
  TF1 *f_Agaus = new TF1("f_Agaus","gaus(0) + pol1(3)",0.7 , 1.1);
  
  f_Agaus->SetParameter(0, 1.);
  f_Agaus->SetParameter(1, 1.);
  f_Agaus->SetParameter(2, 0.01);
  f_Agaus->SetParLimits(2, 0, 1.);
  
  hAoE->Fit("f_Agaus","L0","", 0.7, 1.1);
  
  double A_AoE = f_Agaus->GetParameter(0);
  double mu_AoE = f_Agaus->GetParameter(1);
  double sigma_AoE = f_Agaus->GetParameter(2);
  double sigma_AoE_err = f_Agaus->GetParError(2);
  
  double FWHM_AoE = 2.35*sigma_AoE;
  double FWHM_AoE_err = 2.35*sigma_AoE_err;
  c1->Update();

  ////////////////////////////////////////////////////
  //---------------- PSA analysis --------------------
  ////////////////////////////////////////////////////
  h_AoE2->Add(h_AoE2);
  h_AoE1->Add(h_AoE2);
  h_AoE1->Add(h_AoE3,-1);
  h_AoE->Add(h_AoE1,-1);
  h_AoE->Scale(1./h_AoE->GetBinContent(h_AoE->GetMaximumBin()));
  
  double totIntegral = h_AoE->Integral(1, h_AoE->GetNbinsX());
  cout << "Total integral = " << totIntegral << endl;
  
  double SurvFrac_d = 1.;
  double cutIntegral;
  int bin_lowCut = 1;
  while ( SurvFrac_d >= 0.91 ){
    cutIntegral = h_AoE->Integral(bin_lowCut, h_AoE->GetNbinsX());
    SurvFrac_d = (double)cutIntegral/totIntegral;
    bin_lowCut++;
  }
  double AoE_lowCut = h_AoE->GetBinCenter(bin_lowCut);
  cout << "bin_lowCut = " << bin_lowCut << "   " << "AoE_lowcut = " << AoE_lowCut << endl;      
  cout << "CutIntegral = " << cutIntegral << endl;
  cout << "SF_dep = " << SurvFrac_d << endl;
  cout << "sigma_aoe = " << sigma_AoE << endl;
  double AoE_highCut = mu_AoE + 3.*sigma_AoE;
  cout << "AoE_highcut = " << AoE_highCut << endl;
  hAoE->Draw();  
  
  //----------------------------------------------------
  //                Print the result
  //____________________________________________________
    
  TLine *line_l = new TLine(AoE_lowCut,c1->GetUymin(),AoE_lowCut,1.0);
  line_l->SetLineStyle(10);
  line_l->SetLineColor(2);
  line_l->SetLineWidth(2);
  line_l->Draw();
  TLine *line_h = new TLine(AoE_highCut,c1->GetUymin(),AoE_highCut,1.0);
  line_h->SetLineStyle(10);
  line_h->SetLineColor(2);
  line_h->SetLineWidth(2);
  line_h->Draw();
  TLegend *legend = new TLegend(0.1,0.66,0.52,0.9,"","NDC");
  legend->AddEntry(hAoE,Form("FWHM_{A/E} = (%4.2f #pm %4.2f)%%",FWHM_AoE*100,FWHM_AoE_err*100),"l");
  legend->AddEntry(line_l,Form("A/E_{lowcut} = %5.3f ", AoE_lowCut),"l");
  legend->AddEntry(line_h,Form("A/E_{highcut} = %5.3f ", AoE_highCut),"l");
  legend->Draw();
  c1->Update();
  c1->Print(Form("%s/chn%d_AoE_onlyDEP.pdf",resdir,chn));
  
  
  //-----------Draw the spectrum before and after PSA______________
  TCanvas *c2 = new TCanvas("c2","PSA analysis");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  hene->Draw();
  
  TH1D *hene_psa = new TH1D("hene_psa","; Energy [keV]; Counts", bin_in, mine, maxe);
  for ( int i = 0; i < events; i++){
    double AoEnorm = AoEs.at(i)/mean_aoe;
    if ( AoEnorm <= AoE_highCut && AoEnorm >= AoE_lowCut ) hene_psa->Fill(energies.at(i));
  }
  hene_psa->Draw("same");
  c2->Update();
  
  //------------- fit DEP region after PSD ----------------
  f_2gaus->SetParameter(0, A_dep*0.9);
  f_2gaus->SetParameter(1, mu_dep);
  f_2gaus->SetParameter(2, sigma_dep);
  f_2gaus->SetParLimits(2, 0, 100);
  f_2gaus->SetParameter(3, A_fep*0.15);
  f_2gaus->SetParameter(4, mu_fep);
  f_2gaus->SetParameter(5, sigma_fep);
  f_2gaus->SetParLimits(5, 0, 100);
  f_2gaus->SetParameter(6, bkg*0.5);
  
  gStyle->SetOptFit(0);
  hene_psa->Fit("f_2gaus", "L+", "sames", min, max);
  
  double bkg_psa = f_2gaus->GetParameter(6);
  
  //------------ calculation of Area of DEP,FEP -----------
  //----- DEP Area ------
  
  int min_bin = hene->GetXaxis()->FindBin(mu_dep - n_sigma*sigma_dep);
  int max_bin = hene->GetXaxis()->FindBin(mu_dep + n_sigma*sigma_dep);
  int min_bin_1 = hene->GetXaxis()->FindBin(mu_dep - 2*n_sigma*sigma_dep);
  int max_bin_1 = hene->GetXaxis()->FindBin(mu_dep + 2*n_sigma*sigma_dep);
  int max_bin_2 =  hene->GetXaxis()->FindBin(mu_dep + 3*n_sigma*sigma_dep);
  int delta = max_bin - min_bin;
  
  int count = 0;
  int sumCount = 0;  
  int count_psa = 0;
  int sumCount_psa = 0;
  for (Int_t j = min_bin; j <= max_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;
    count_psa = hene_psa->GetBinContent(j);
    sumCount_psa += count_psa;
  }
  double Area_d_tot = sumCount;
  double Area_d_psa_tot = sumCount_psa;
  cout << "AreaDEP_tot = " << Area_d_tot << endl;
  cout << "Area_DEP_tot_psa = " << Area_d_psa_tot << endl;
  
  //---background
  
  count = 0, sumCount = 0, count_psa = 0, sumCount_psa = 0;
  for (Int_t j =  min_bin_1; j <= min_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;
    count_psa = hene_psa->GetBinContent(j);
    sumCount_psa += count_psa;
  }
  double b1_I = sumCount;
  double b2_I = sumCount_psa;
  count = 0, sumCount = 0, count_psa = 0, sumCount_psa = 0;
  for (Int_t j =  max_bin; j <= max_bin_1; j++){
    count = hene->GetBinContent(j);
    sumCount += count;
    count_psa = hene_psa->GetBinContent(j);
    sumCount_psa += count_psa;
  }
  double b1_II = sumCount;
  double b2_II = sumCount_psa;
  count = 0, sumCount = 0, count_psa = 0, sumCount_psa = 0;
  for (Int_t j =  max_bin_1; j <= max_bin_2; j++){
    count = hene->GetBinContent(j);
    sumCount += count;
    count_psa = hene_psa->GetBinContent(j);
    sumCount_psa += count_psa;
  }    
  double b1_III = sumCount;
  double b2_III = sumCount_psa;
  
  double b1_d = b1_I + 2*b1_II - b1_III;
  double Area_d = Area_d_tot - b1_d;
  cout << "AreaDEP = " << Area_d << endl;
  double b2_d = b2_I + 2*b2_II - b2_III;
  double Area_d_psa = Area_d_psa_tot - b2_d;//subtraction of background
  cout << "Area_DEP_psa = " << Area_d_psa << endl;

  //////////////////////
  //----- FEP area -----
  //////////////////////
  
  min_bin = hene->GetXaxis()->FindBin(mu_fep - n_sigma*sigma_fep);
  max_bin = hene->GetXaxis()->FindBin(mu_fep + n_sigma*sigma_fep);
  min_bin_1 = hene->GetXaxis()->FindBin(mu_fep - 2*n_sigma*sigma_fep);
  max_bin_1 = hene->GetXaxis()->FindBin(mu_fep + 2*n_sigma*sigma_fep);
  delta = max_bin - min_bin;
  
  count = 0, sumCount = 0, count_psa = 0, sumCount_psa = 0;
  for (Int_t j = min_bin; j <= max_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;
    count_psa = hene_psa->GetBinContent(j);
    sumCount_psa += count_psa;
  }
  double Area_f_tot = sumCount;
  double Area_f_psa_tot = sumCount_psa;
  cout << "Area_FEP_tot = " << Area_f_tot << endl;
  cout << "Area_FEP_tot_psa = " << Area_f_psa_tot << endl;
  
  //---background
  
  count = 0, sumCount = 0, count_psa = 0, sumCount_psa = 0;
  for (Int_t j =  min_bin_1; j <= min_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;
    count_psa = hene_psa->GetBinContent(j);
    sumCount_psa += count_psa;
  }
  for (Int_t j =  max_bin; j <= max_bin_1; j++){
    count = hene->GetBinContent(j);
    sumCount += count;
    count_psa = hene_psa->GetBinContent(j);
    sumCount_psa += count_psa;
  }
  double b1_f = sumCount;
  double Area_f = Area_f_tot - b1_f;
  cout << "Area_FEP = " << Area_f << endl;
  double b2_f = sumCount_psa;
  double Area_f_psa = Area_f_psa_tot - b2_f;
  cout << "Area_FEP_psa = " << Area_f_psa << endl;
  
  //------- Acceptance and error calculation --------
  
  double Acc_d = Area_d_psa/Area_d;
  double Acc_f = Area_f_psa/Area_f;
  double Acc_bkg = bkg_psa/bkg;
  
  double Acc_d_e = pow((Area_d_psa_tot+b2_d)*Area_d*Area_d + Area_d_psa*Area_d_psa*(Area_d_tot+b1_d),0.5)/(Area_d*Area_d);
  double Acc_f_e = pow((Area_f_psa_tot+b2_f)*Area_f*Area_f + Area_f_psa*Area_f_psa*(Area_f_tot+b1_f),0.5)/(Area_f*Area_f);
    
  //----- print the results of PSD -----

  hene->SetFillColor(16);//color of the spectrum before PSD
  hene->SetLineColor(16);
  hene_psa->SetFillColor(4);//color of the spectrum after PSD 
  hene_psa->SetLineColor(4);

  TLegend *leg = new TLegend(0.57,0.76,0.99,0.99,"","NDC");
  leg->AddEntry(hene,"Before PSD","f");
  leg->AddEntry(hene_psa,"After PSD","f");
  leg->AddEntry("0",Form("SF_{DEP} = (%4.2f #pm %4.2f)%% ",Acc_d*100,Acc_d_e*100),"0");
  leg->AddEntry("0",Form("SF_{FEP} = (%4.2f #pm %4.2f)%% ",Acc_f*100,Acc_f_e*100),"0");
  leg->Draw();
  c2->Update();
  c2->Print(Form("%s/chn%d_PSA_DEP_FEP.pdf",resdir,chn));
  
   
  ////////////////////////////////////////////////////////////////
  //------------------- AoE correction ---------------------------
  ////////////////////////////////////////////////////////////////
    
  //vector<double> AoE_corr(TH1D hist, int chn, double min_aoe_dep,  double calib, int shapingMethod, int shap){
  cout << "Correction of AoE..." << endl;
    
  //------ energy to analize -----
  double ene[20], ene_err[20], mu_AoE_[20], mu_AoE_err_[20];
  ene[0] = 900.;
  ene[1] = 1000.;
  ene[2] = 1150.;
  ene[3] = 1250.;
  ene[4] = 1375.;
  ene[5] = 1425.;
  ene[6] = 1550.;
  ene[7] = 1720.;
  ene[8] = 1845.;
  ene[9] = 1900.;
  ene[10] = 1955.;
  ene[11] = 2010.;
  //ene[8] = 2065.;
  ene[12] = 2155.;
  ene[13] = 2245.;
  ene[14] = 2300.;
  ene[15] = 2355.5;
  int nPeaks = 16;
    
  TCanvas *c3 = new TCanvas("c3","AoE-correction");
  c3->Divide(4,4);

  for(int j = 0; j < nPeaks; j++){
    //--------- find parameters from spectrum -------
    ene_err[j]=0;
    double peak = ene[j];
    double min = peak - 20;
    double max = peak + 20;
    int bin = ( max - min ) * 2 ;
    
    c3->cd(j+1);
    gStyle->SetPalette(1);
    
    TH1D *h_AoE = new TH1D("h_AoE",Form("AoE - region %g;A/E;Counts",chn,ene[j]), 200, 0.7, 1.1);
    for ( int i = 0; i < events; i++){
      double AoEnorm = AoEs.at(i)/mean_aoe;
      if ( energies.at(i) <= max &&  energies.at(i) >= min ) h_AoE->Fill(AoEnorm);
    }
    h_AoE->Draw();
    c3->Update();
    
    //---------------- fit the A/E ------------
    
    double fMin = 0;
    double fMax = 3;
    
    TF1 *FitAoE = new TF1( "FitAoE","[0]*TMath::Exp( -TMath::Power( (x-[1])/[2],2. )/2. ) +[3]*((TMath::Exp([4]*(x-[5])))+[6])/((TMath::Exp((x-[5])/[7]))+[5])",fMin, fMax );
    
    double A_AoE = h_AoE->GetMaximum();
    double mu_AoE = h_AoE->GetBinCenter( h_AoE->GetMaximumBin());
    double sigma_AoE = mu_AoE*0.01;
    
    FitAoE->SetParameter(0, A_AoE*0.9);
    FitAoE->SetParLimits(0, A_AoE*0.6, A_AoE*1.2);
    FitAoE->SetParameter(1, mu_AoE);
    FitAoE->SetParLimits(1, 0.98, 1.02);
    FitAoE->SetParameter(2, 0.01);
    FitAoE->SetParLimits(2, 0.005, 0.03);
    FitAoE->SetParameter( 3, A_AoE/3. );
    FitAoE->SetParameter( 4, 0.1 );
    //FitAoE->SetParLimits( 4, 0, 1);
    FitAoE->SetParameter( 5, 0.98 );
    //FitAoE->SetParLimits( 5, 0.7, 1. );
    FitAoE->SetParameter( 6, -0.8 );
    FitAoE->SetParLimits( 6, -1, 0 );
    FitAoE->SetParameter( 7, 0.02);
    FitAoE->SetParLimits( 7, 0.005, 0.1 );
    
    gStyle->SetOptFit(0);
    h_AoE->Fit("FitAoE","L","", 0.7, 1.1);
    TF1 *gaussian = new TF1("gaussian","[0]*TMath::Exp( -TMath::Power( (x-[1])/[2],2. )/2. )",0.7,1.1);
    gaussian->SetParameters(FitAoE->GetParameter(0),FitAoE->GetParameter(1),FitAoE->GetParameter(2));
    
    TF1 *step = new TF1("step", "[0]*((TMath::Exp([1]*(x-[2])))+[3])/((TMath::Exp((x-[2])/[4]))+[2])",fMin, fMax );
    step->SetParameters(FitAoE->GetParameter(3), FitAoE->GetParameter(4), FitAoE->GetParameter(5), FitAoE->GetParameter(6), FitAoE->GetParameter(7));
    
    gaussian->SetLineColor(kAzure);
    step->SetLineColor(kSpring+2);
    gaussian->SetLineStyle(2);
    step->SetLineStyle(2);
    gaussian->Draw("same");
    step->Draw("same");
    c3->Update();
    
    mu_AoE_[j] = FitAoE->GetParameter(1);
    mu_AoE_err_[j] = FitAoE->GetParError(1);
  }
  c3->Print(Form("%s/chn%d_AoE_fits.pdf",resdir,chn));
  
  //----- AoE versus E ----
    
  TCanvas *c_3 = new TCanvas("c_3","AoE vs E");
  TGraph *g_e = new TGraphErrors (nPeaks,ene,mu_AoE_,ene_err,mu_AoE_err_);
    
  g_e->SetTitle(Form("chn%d  #mu_{A/E} vs E ",chn));
  g_e->GetXaxis()->SetTitle("Energy [keV]");
  g_e->GetYaxis()->SetTitle("");
  g_e->SetMarkerColor(4);
  g_e->SetLineColor(4);
  g_e->SetMarkerStyle(21);
  g_e->Draw("AP");
    
  TF1 *pol1 = new TF1("pol1","pol1",600,2400);
  gStyle->SetOptFit(111);

  pol1->SetParameter(0, 1);
  pol1->SetParameter(1, pow(10,-6.));
  pol1->SetLineColor(kRed+1);

  g_e->Fit("pol1","","",600,2400);
  c_3->Update();
    
  double par1 = pol1->GetParameter(0);
  double par2 = pol1->GetParameter(1);
  double errpar1 = pol1->GetParError(0);
  double errpar2 = pol1->GetParError(1);

  TF1 *fit_plus = new TF1("fit_plus",Form("%g + %g*x + pow((%g*%g + x*x*%g*%g),0.5)",par1,par2,errpar1,errpar1,errpar2,errpar2),600,2400);
  fit_plus->SetLineColor(kRed-7);
  fit_plus->SetLineStyle(2);
  fit_plus->Draw("same");

  TF1 *fit_minus = new TF1("fit_minus",Form("%g + %g*x - pow((%g*%g + x*x*%g*%g),0.5)",par1,par2,errpar1,errpar1,errpar2,errpar2),600,2400);
  fit_minus->SetLineColor(kRed-7);
  fit_minus->SetLineStyle(2);
  fit_minus->Draw("same");
  c_3->Update();
  c_3->Print(Form("%s/chn%d_AoE_vs_E.pdf",resdir,chn));

  vector<double> b;
  b.push_back(par1);
  b.push_back(par2);
  b.push_back(errpar1);
  b.push_back(errpar2);
  // b.push_back(min_AoE);
  //return(b);
  for (int i = 0; i < 4; i++)
     cout << b.at(i) << endl;

  
  ////////////////// ////////////////////////////////////////////////
  //------------------Survival fraction of SEP----------------------
  //////////////////////////////////////////////////////////////////

  TCanvas *c_4 = new TCanvas("c_4","energy");
  gPad->SetLogy();

  double peak = 2104.;
  double min_AoE_corr = AoE_lowCut + b.at(1)*(peak-1592.5);
  double max_AoE_corr = AoE_highCut + b.at(1)*(peak-1592.5);
  cout << "Calculation of survival fraction for " << peak << endl;

  //------- peak region ----------                                   
  min = (peak - 45);
  max = (peak + 45);
  hene->Draw();
  hene->GetXaxis()->SetRangeUser(min,max);
  c_4->Update();

  //----- fit of the energy peaks with a gaussian --------
  TF1 *f_gaus = new TF1("f_gaus","gaus(0)+ (1. - 0.5*TMath::Erfc((x - [1])/[2]/TMath::Sqrt(2.0))) + pol0(3)",min,max);

  //--------- initialization of parameters of the fit ---------
  double A = hene->GetMaximum();
  double mu = hene->GetBinCenter(hene->GetMaximumBin());
  f_gaus->SetParameter(0,A);
  f_gaus->SetParameter(1,mu);
  f_gaus->SetParameter(2,mu*0.0005);
  f_gaus->SetParLimits(2,0,100);

  hene->GetXaxis()->SetRangeUser(2060,2090);
  max_bkg = hene->GetMaximum();
  min_bkg = hene->GetMinimum();
  double bkg_2104 = (max_bkg+min_bkg)/2.;
  f_gaus->SetParameter(3,bkg_2104);
  
  cout << "Initialization of bkg_2104 =  " << bkg_2104 << endl;
  hene->GetXaxis()->SetRangeUser(min,max);
  gStyle->SetOptFit(0);
  hene->Fit("f_gaus","L","",min,max);
  c_4->Update();

  double A_ = f_gaus->GetParameter(0);
  double mu_ = f_gaus->GetParameter(1);
  double sigma_ = f_gaus->GetParameter(2);
  bkg = f_gaus->GetParameter(3);

  //----draw the A/E----    
  TH1D *h_sep = new TH1D("h_sep","", 100, 0.7, 1.1);
  for ( int i = 0; i < events; i++){
    double AoEnorm = AoEs.at(i)/mean_aoe;
    if ( energies.at(i) <= mu_ + n_sigma*sigma_ &&  energies.at(i) >= mu_ - n_sigma*sigma_ ) h_sep->Fill(AoEnorm);
  }
  h_sep->Scale(1./h_sep->GetBinContent(h_sep->GetMaximumBin()));
  
  //------------ calculation of Area of peak ---------
  min_bin = hene->GetXaxis()->FindBin(mu_ - n_sigma * sigma_);
  max_bin = hene->GetXaxis()->FindBin(mu_ + n_sigma * sigma_);
  min_bin_1 = hene->GetXaxis()->FindBin(mu_ - 2*n_sigma * sigma_);
  max_bin_1 = hene->GetXaxis()->FindBin(mu_ + 2*n_sigma * sigma_);
  delta = max_bin - min_bin;
  
  count = 0, sumCount = 0;
  for (Int_t j = min_bin; j <= max_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;}
  double Area_tot_sep = sumCount;//total area                                     
  cout << "Area_tot_sep = " << Area_tot_sep << endl;

  //---background                                                          
  count = 0, sumCount = 0;
  for (Int_t j =  min_bin_1; j <= min_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;}
  for (Int_t j =  max_bin; j <= max_bin_1; j++){
    count = hene->GetBinContent(j);
    sumCount += count;}
  double b1_sep = sumCount;
  double Area_sep = Area_tot_sep - b1_sep;//subtraction of background                                        
  cout << "Area_sep = " << Area_sep << endl;

  //---------------- PSA analysis -------------------- 
  hene_psa = new TH1D("hene_psa","; Energy [keV]; Counts",bin_in, mine, maxe);
  c_4->cd();
  for ( int i = 0; i < events; i++){
    double AoEnorm = AoEs.at(i)/mean_aoe;
    if ( AoEnorm <= max_AoE_corr && AoEnorm >= min_AoE_corr ) hene_psa->Fill(energies.at(i));
  }
  hene_psa->Draw("same");
  hene_psa->GetXaxis()->SetRangeUser(min,max);
  c_4->Update();

  //------------- fit SEP region after PSD ----------------                         
  f_gaus->SetParameter(0, A_*0.3);
  f_gaus->SetParameter(1, mu_);
  f_gaus->SetParameter(2, sigma_);
  f_gaus->SetParLimits(2, 0, 100);
  f_gaus->SetParameter(3, bkg*0.5);

  gStyle->SetOptFit(0);
  hene_psa->Fit("f_gaus","L","sames",min, max);

  //-------------- get the results from fit --------------               
  double mu_psa = f_gaus->GetParameter(1);
  double sigma_psa = f_gaus->GetParameter(2);

  //------------ calculation of Area of peak -----------
  count = 0, sumCount = 0;
  for (Int_t j = min_bin; j < max_bin + 1; j++){
     count = hene_psa->GetBinContent(j);
     sumCount += count;}
  double Area_psa_tot_sep = sumCount;//total area
  cout << "Area_tot_psa_sep = " << Area_psa_tot_sep << endl;

  //---background
  count = 0, sumCount = 0;
  for (Int_t j =  min_bin_1; j <= min_bin; j++){
    count = hene_psa->GetBinContent(j);
    sumCount += count;}
  for (Int_t j =  max_bin; j <= max_bin_1; j++){
    count = hene_psa->GetBinContent(j);
    sumCount += count;}
  double b2_sep = sumCount;
  double Area_psa_sep = Area_psa_tot_sep - b2_sep;//subtraction of background
  cout << "Area_psa_sep = " << Area_psa_sep << endl;

  //------- Accentance calculation --------
  double Acc_sep = 1.*Area_psa_sep/Area_sep;
  double Acc_e_sep = pow((Area_psa_tot_sep+b2_sep)*Area_sep*Area_sep + Area_psa_sep*Area_psa_sep*(Area_tot_sep+b1_sep),0.5)/(Area_sep*Area_sep);

  //--------------------------------------------
  //          print the results of PSD 
  //____________________________________________                  
  c_4->cd();
  hene->SetFillColor(16);//color of the spectrum befor psd
  hene->SetLineColor(16);
  hene_psa->SetFillColor(4);//color of the spectrum after psd
  hene_psa->SetLineColor(4);

  TLegend *legen_sep = new TLegend(0.52,0.8,0.99,0.99,"","NDC");
  legen_sep->AddEntry(hene,"Before PSD","f");
  legen_sep->AddEntry(hene_psa,"After PSD","f");
  legen_sep->AddEntry("0",Form("SF_{%g} = (%4.2f #pm %4.2f)%%",peak, Acc_sep*100,Acc_e_sep*100),"0");
  legen_sep->Draw();
  c_4->Update();
  c_4->Print(Form("%s/chn%d_PSA_peak%g.pdf",resdir,chn,peak));


  //////////////////////////////////////////////////////
  //-----------Survival fraction of FEP_2615-----------
  /////////////////////////////////////////////////////

  TCanvas *c_5 = new TCanvas("c_5","energy");
  gPad->SetLogy();

  TH1D *h_2614 = new TH1D("h_2614","", 100, 0.7, 1.1);;

  peak = 2614.5;
  min_AoE_corr = AoE_lowCut + b.at(1)*(peak-1592.5);
  max_AoE_corr = AoE_highCut + b.at(1)*(peak-1592.5); 
  cout << "Calculation of survival fraction for " << peak << endl;

  //------- peak region ----------  
  min = (peak - 45);
  max = (peak + 45);
  hene->Draw();
  hene->GetXaxis()->SetRangeUser(min,max);
  c_5->Update();

  //----- fit of the energy peaks with a gaussian and an error function --------
  TF1 *fit_fep = new TF1( "fit_fep","gaus(0)+0.5*TMath::Erfc((x - [1])/[2]/TMath::Sqrt(2.0)) + pol0(3)",min,max);

  //--------- initialization of parameters of the fit ----------
  A = hene->GetMaximum();
  mu = hene->GetBinCenter(hene->GetMaximumBin());
  fit_fep->SetParameter(0,A);
  fit_fep->SetParameter(1,mu);
  fit_fep->SetParameter(2,mu*0.0005);
  fit_fep->SetParLimits(2,0,100);
  fit_fep->SetParameter(3, 1.);  

  hene->GetXaxis()->SetRangeUser(2570,2600);
  max_bkg = hene->GetMaximum();
  min_bkg = hene->GetMinimum();
  double bkg_2615 = (max_bkg+min_bkg)/2.;
  fit_fep->SetParameter(3,bkg_2104);

  cout << "Initialization of bkg_2615 = " << bkg_2615 << endl;
  hene->GetXaxis()->SetRangeUser(min,max);
  gStyle->SetOptFit(0);
  hene->Fit("fit_fep","L","",min,max);
  c_5->Update();

  A_ = fit_fep->GetParameter(0);
  mu_ = fit_fep->GetParameter(1);
  sigma_ = fit_fep->GetParameter(2);
  bkg = fit_fep->GetParameter(3);

  //---- A/E----
  for ( int i = 0; i < events; i++){
    double AoEnorm = AoEs.at(i)/mean_aoe;
    if ( energies.at(i) <= mu_ + n_sigma*sigma_ &&  energies.at(i) >= mu_ - n_sigma*sigma_ ) h_2614->Fill(AoEnorm);
  }
  h_2614->Scale(1./h_2614->GetBinContent(h_2614->GetMaximumBin()));

  //------------ calculation of Area of peak ---------
  min_bin = hene->GetXaxis()->FindBin(mu_ - n_sigma * sigma_);
  max_bin = hene->GetXaxis()->FindBin(mu_ + n_sigma * sigma_);
  min_bin_1 = hene->GetXaxis()->FindBin(mu_ - 2*n_sigma * sigma_);
  max_bin_1 = hene->GetXaxis()->FindBin(mu_ + 2*n_sigma * sigma_);
  delta = max_bin - min_bin;

  count = 0, sumCount = 0;
  for (Int_t j = min_bin; j <= max_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;}
  double Area_tot_fep = sumCount;//total area                                                                            
  cout << "Area_tot_fep = " << Area_tot_fep << endl;

  //---background 
  count = 0, sumCount = 0;
  for (Int_t j =  min_bin_1; j <= min_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;}
  for (Int_t j =  max_bin; j <= max_bin_1; j++){
    count = hene->GetBinContent(j);
    sumCount += count;}
  double b1_fep = sumCount;
  double Area_fep = Area_tot_fep - b1_fep;//subtraction of background                                                                                  
  cout << "Area_fep = " << Area_fep << endl;

  //////////////////////////////////////////////////// 
  //---------------- PSA analysis --------------------

  hene_psa = new TH1D("hene_psa","; Energy [keV]; Counts",bin_in, mine, maxe);
  c_5->cd();
  for ( int i = 0; i < events; i++){
    double AoEnorm = AoEs.at(i)/mean_aoe;
    if ( AoEnorm <= max_AoE_corr && AoEnorm >= min_AoE_corr ) hene_psa->Fill(energies.at(i));
  }
  hene_psa->Draw("same");
  hene_psa->GetXaxis()->SetRangeUser(min,max);
  c_5->Update();

  //------------- fit DEP region after PSD ---------------
  fit_fep->SetParameter(0, A_*0.3);
  fit_fep->SetParameter(1, mu_);
  fit_fep->SetParameter(2, sigma_);
  fit_fep->SetParLimits(2, 0, 100);
  fit_fep->SetParameter(3, bkg*0.5);

  gStyle->SetOptFit(0);
  hene_psa->Fit("fit_fep","L","sames",min, max);

  //-------------- get the results from fit --------------
  mu_psa = fit_fep->GetParameter(1);
  sigma_psa = fit_fep->GetParameter(2);

  //------------ calculation of Area of peak -----------
  count = 0, sumCount = 0;
  for (Int_t j = min_bin; j < max_bin + 1; j++){
    count = hene_psa->GetBinContent(j);
    sumCount += count;}
  double Area_psa_tot_fep = sumCount;//total area
  cout << "Area_tot_psa_fep = " << Area_psa_tot_fep << endl;

  //---background
  count = 0, sumCount = 0;
  for (Int_t j =  min_bin_1; j <= min_bin; j++){
    count = hene_psa->GetBinContent(j);
    sumCount += count;}
  for (Int_t j =  max_bin; j <= max_bin_1; j++){
    count = hene_psa->GetBinContent(j);
    sumCount += count;}
  double b2_fep = sumCount;
  double Area_psa_fep = Area_psa_tot_fep - b2_fep;//subtraction of background         
  cout << "Area_psa_fep = " << Area_psa_fep << endl;

  //------- Accentance calculation --------

  double Acc_fep = 1.*Area_psa_fep/Area_fep;
  double Acc_e_fep = pow((Area_psa_tot_fep+b2_fep)*Area_fep*Area_fep + Area_psa_fep*Area_psa_fep*(Area_tot_fep+b1_fep),0.5)/(Area_fep*Area_fep);

  //--------------------------------------------
  //          print the results of PSD      
  //____________________________________________ 

  c_5->cd();
  hene->SetFillColor(16);//color of the spectrum befor psd
  hene->SetLineColor(16);
  hene_psa->SetFillColor(4);//color of the spectrum after psd
  hene_psa->SetLineColor(4);

  TLegend *legen_fep = new TLegend(0.52,0.8,0.99,0.99,"","NDC");
  legen_fep->AddEntry(hene,"Before PSD","f");
  legen_fep->AddEntry(hene_psa,"After PSD","f");
  legen_fep->AddEntry("0",Form("SF_{%g} = (%4.2f #pm %4.2f)%%",peak, Acc_fep*100,Acc_e_fep*100),"0");
  legen_fep->Draw();
  c_5->Update();
  c_5->Print(Form("%s/chn%d_PSA_peak%g.pdf",resdir,chn,peak));


  //////////////////////////////////////////////////////
  //----------- survival fraction background -----------
  //////////////////////////////////////////////////////
    
  double qbb = 2039.;
  cout << "Calculation of survival fraction background at " << qbb << endl;
  min_AoE_corr = AoE_lowCut + b.at(1)*(qbb-1592.5);
  max_AoE_corr = AoE_highCut + b.at(1)*(qbb-1592.5); 

  //------- qbb region ----------
  min = (qbb - 35);
  max = (qbb + 35);
  c_4->cd();
  hene->Draw();
  hene->GetXaxis()->SetRangeUser(min,max);

  //----------- A/E -----------
  TH1D *h_bkg = new TH1D("h_bkg","", 100, 0.7, 1.1);
  for ( int i = 0; i < events; i++){
    double AoEnorm = AoEs.at(i)/mean_aoe;
    if ( energies.at(i) <= max && energies.at(i) >= min ) h_bkg->Fill(AoEnorm);
  }
  h_bkg->Scale(1./h_bkg->GetBinContent(h_bkg->GetMaximumBin()));
    

  //------------ calculation of Area -----------
  min_bin = hene->GetXaxis()->FindBin(min);
  max_bin = hene->GetXaxis()->FindBin(max);
  sumCount = 0;
  for (Int_t j = min_bin; j <= max_bin; j++){
    count = hene->GetBinContent(j);
    sumCount += count;
  }
  double b2039 = sumCount;//bkg
    
  //---------------- PSA --------------------
  hene_psa = new TH1D("hene_psa","; Energy [keV]; Counts",bin_in,mine,maxe);
  for ( int i = 0; i < events; i++){
    double AoEnorm = AoEs.at(i)/mean_aoe;
    if ( AoEnorm <= max_AoE_corr && AoEnorm >= min_AoE_corr ) hene_psa->Fill(energies.at(i));
  }
  c_4->cd();
  hene->Draw();
  hene->GetXaxis()->SetRangeUser(min,max);
  hene_psa->Draw("same");
  c_4->Update();
 
  //------------ calculation of Area -----------
  sumCount = 0;
  for (Int_t j = min_bin; j <= max_bin; j++){
    count = hene_psa->GetBinContent(j);
    sumCount += count;
  }
  double b2039_psa = sumCount;//total area
  cout << "Area_tot_psa = " << bkg_psa << endl;
    
  //------- Accentance calculation --------
  double SF_bkg = 1.*b2039_psa/b2039;
  double SF_bkg_err = pow(((b2039 + b2039_psa)*b2039_psa)/pow(b2039, 3.), 0.5);
    
  //----- Print the results -----
  
  hene->SetFillColor(16);//color of the spectrum before psd                                                               
  hene->SetLineColor(16);
  hene_psa->SetFillColor(4);//color of the spectrum after psd                                                             
  hene_psa->SetLineColor(4);
 
  TLegend *leg_bkg = new TLegend(0.52,0.8,0.99,0.99,"","NDC");
  leg_bkg->AddEntry(hene,"Before PSD","f");
  leg_bkg->AddEntry(hene_psa,"After PSD","f");
  leg_bkg->AddEntry("0",Form("SF_{ROI} = (%4.2f #pm %4.2f)%%", SF_bkg*100,SF_bkg_err*100),"0");
  leg_bkg->Draw();
  c_4->Update();
  c_4->Print(Form("%s/chn%d_PSA_bkg2039.pdf",resdir,chn));

  //----------------AoE Print------------
  c1->cd();
  gStyle->SetOptStat(0);
  h_AoE->Draw();
  h_AoE->SetLineColor(4);
  h_AoE->SetLineWidth(2);
  h_sep->Draw("same");
  h_sep->SetLineColor(1);
  h_sep->SetLineWidth(2);
  h_2614->Draw("same");
  h_2614->SetLineColor(3);
  h_2614->SetLineWidth(2);
  h_bkg->Draw("same");
  h_bkg->SetLineColor(2);
  h_bkg->SetLineWidth(2);

  c1->Update();
  TLegend *leg_p = new TLegend(0.1,0.7,0.3,0.9,"","NDC");
  leg_p->AddEntry(h_AoE,"DEP","l");
  leg_p->AddEntry(h_sep,"SEP","l");
  leg_p->AddEntry(h_2614,"FEP_{2615}","l");
  leg_p->AddEntry(h_bkg,"ROI","l");
  leg_p->Draw();
  c1->Update();
  c1->Print(Form("%s/chn%d_AoE.pdf",resdir,chn));
    
  //-------- Save results in file --------
  ofstream file;
  char *out_file = Form("%s/PSA_results.txt",resdir);
  char *out_results = Form("%d %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n",chn, FWHM_AoE*100, FWHM_AoE_err*100, Acc_d*100, Acc_d_e*100, Acc_f*100, Acc_f_e*100,Acc_sep*100, Acc_e_sep*100,Acc_fep*100, Acc_e_fep*100, SF_bkg*100, SF_bkg_err*100);
    
  file.open(out_file,std::ios_base::app|std::ios_base::out);
  file << out_results;
  file.close();
  //myapp->Run();
  return 0;
}
