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
#include "TF1.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include <TApplication.h>

using namespace std;

int main( int argc, char* argv[]){
  //TApplication * myapp = new TApplication("myapp",0,0);
  char *resdir = argv[1];
  double cycle = 6.00;
  if (argc>2) cycle = atof(argv[2]);
  int run = 98;
  if (argc>3) run = atoi(argv[3]);
  int weight = 10;
  if (argc>4) weight = atoi(argv[4]);
  
  cout << "Read files from run " << run << " with weight " << weight << endl;
  
  //char *filelist = Form("/nfs/gerda5/gerda-data/blind/v0%3.2f/meta/data-sets/cal/run0095-run0114-cal-analysis.txt",cycle);
  char *filelist = Form("/nfs/gerda5/gerda-data/blind/v0%3.2f/meta/data-sets/cal/run00%d-cal-analysis.txt", cycle, run);
  cout << "Read events from " << filelist << endl;  
  ifstream list(filelist);
  
  int totEntries = 0;
  
  const int nChn = 41;
  
  char *filename = Form("%s/run00%d-cal_dplms_%d.tier.root", resdir, run, weight);
  cout << "Writing on file " << filename << endl;
  
  TFile *newfile = new TFile(filename,"RECREATE");
  //TFile *newfile = new TFile(Form("%s/run0095-run0114-cal.tier.root",resdir),"RECREATE");
  
  TTree *tier = new TTree("tier","Tree with energy and AoE");
  ULong64_t tempo;
  tier->Branch("tempo", &tempo);
  std::vector<double>* energy = NULL;
  tier->Branch("energy", &energy);
  std::vector<double>* AoE = NULL;
  tier->Branch("AoE", &AoE);
  std::vector<double>* AoE_dplms = NULL;
  tier->Branch("AoE_dplms", &AoE_dplms);
  
  
  int count[nChn] = {0};
  clock_t begin, now, time;
  begin = clock();//initial time
  while (1){
    //for (int conteggi = 0; conteggi < 2; conteggi++){
    char line[150];
    list.getline(line, 128);
    if (!list.good()) break;
    
    char trun[7];
    memcpy(trun, line + 6,7);
    cout << "Run number = " << trun << endl;
    
    //char *tier2file = Form("/nfs/gerda5/gerda-data/blind/v0%3.2f/gen/tier2/ged/cal/run00%d/%s-ged-tier2.root", cycle, run, line);
    //char *tier3file = Form("/nfs/gerda5/gerda-data/blind/v0%3.2f/gen/tier3/all/cal/run00%d/%s-all-tier3.root", cycle, run, line);
    
    //char *tier2file = Form("/nfs/gerda5/gerda-data/blind/v0%3.2f/gen/tier2/ged/cal/%s/%s-ged-tier2.root", cycle, trun, line);
    
    char *tier2file = Form("/nfs/gerda6/users/dandrea/Analysis/createPSDfilter/%s/tier2/weight_%d/%s-ged-tier2.root", trun, weight, line);
    char *tier3file = Form("/nfs/gerda5/gerda-data/blind/v0%3.2f/gen/tier3/all/cal/%s/%s-all-tier3.root", cycle, trun, line);
    cout << "Using file: " << tier2file << endl;
    cout << "Using file: " << tier3file << endl;
    
    TFile *tier2_base = TFile::Open(tier2file);
    TTree *tier2 = (TTree*) tier2_base->Get("tier2");
    
    TFile *tier3_base = TFile::Open(tier3file);
    TTree *tier3 = (TTree*) tier3_base->Get("tier3");
    
    //tier3
    std::vector<double>* energycal = NULL;    
    tier3->SetBranchAddress("energy", &energycal);
    //std::vector<double>* AsuE = NULL;    
    //tier3->SetBranchAddress("AoE", &AsuE);    
    std::vector<int>* failedFlag = NULL;
    tier3->SetBranchAddress("failedFlag", &failedFlag);
    int isTP = -1;    
    tier3->SetBranchAddress("isTP", &isTP);
    Int_t eventChannelNumber;
    tier3->SetBranchAddress("eventChannelNumber", &eventChannelNumber);
    ULong64_t timestamp;
    UInt_t decimalTimestamp;
    tier3->SetBranchAddress ("timestamp",  &timestamp);
    tier3->SetBranchAddress ("decimalTimestamp",  &decimalTimestamp);

    //tier2
    std::vector<double> *A = new vector<double>;
    tier2->SetBranchAddress ("GEMDCurrentPSA_A",  &A);
    std::vector<double> *A_dplms = new vector<double>;
    //tier2->SetBranchAddress ("GEMDGenericShaping_A_energy",  &A_dplms);
    tier2->SetBranchAddress ("GEMDMinMaxFinder_A_maxAmp",  &A_dplms);
    std::vector<double> *energyGauss = new vector<double>;
    tier2->SetBranchAddress ("GEMDEnergyGauss_energy",  &energyGauss);
    //std::vector<double> *energyZAC = new vector<double>;
    //tier2->SetBranchAddress ("GEMDZACShaping_energy",  &energyZAC);
    
    const int nEntries = tier2->GetEntries();
    totEntries += nEntries;
    cout << "Number of events = " << nEntries << endl;
    
    for ( int iEntry = 0; iEntry < nEntries; iEntry++ ) {
      if ((iEntry+1) % 10000 == 0) cout << "event n." << iEntry+1 << endl;
      tier2->GetEntry(iEntry);
      tier3->GetEntry(iEntry);
      energy->clear();
      AoE->clear();
      AoE_dplms->clear();
      tempo = timestamp;
      for (int i = 0; i < eventChannelNumber; i++) {
	double ene = energycal->at(i);
	//if ( !isTP && !failedFlag->at(i) && i==chn && ene>1000 ){
	if ( !isTP && !failedFlag->at(i) && ene>50 && ene<3000 ){
	  count[i]++;
	  //double AoE_official = AoE->at(i);
	  //double eneZAC = energyZAC->at(i);
	  double eneGauss = energyGauss->at(i);
	  double current = A->at(i);
	  double current_dplms = A_dplms->at(i);
	  energy->push_back(ene);
	  AoE->push_back(current/eneGauss);
	  AoE_dplms->push_back(current_dplms/eneGauss);
	  
	  if( count[i] < 1 ) cout << i << " energy = " << ene << " Gauss = " << eneGauss << " A = " << current << " A-dplms = " << current_dplms << endl;
	}
	else {
	  energy->push_back(0);
	  AoE->push_back(0);
	  AoE_dplms->push_back(0);
	}
      }
      tier->Fill();
    }//loop sulle entries
    tier2_base->Close();
    tier3_base->Close();
  }
  newfile->Write();
  now = clock();
  time = now - begin;
  cout << "The time to files is " << time/1.e6 << " s" << endl;
  
  for (int i = 0; i < nChn; i++)
    cout << "chn " << i << ", Total events = " << count[i] << endl;
  
  //myapp->Run();
  return 0;
}
