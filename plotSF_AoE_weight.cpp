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
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TTree.h"
#include "TLatex.h"
#include "TText.h"
#include "TCut.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TPaveLabel.h"

using namespace std;
int main( int argc, char* argv[]){
  //TApplication * myapp = new TApplication("myapp",0,0);                                           
  char *name = argv[1];
  char *resdir = argv[2];
  int chn = atoi(argv[3]);
 
  ifstream file;
  file.open(name);
  
  double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
  double x11, x12, x13, x14, x15, x16;
  double AoE[1000];
  double AoE_err[1000];
  double SF_dep[1000];
  double SF_dep_err[1000];
  double SF_fep[1000];
  double SF_fep_err[1000];
  double SF_sep[1000];
  double SF_sep_err[1000];
  double SF_2614[1000];
  double SF_2614_err[1000];
  double SF_bkg[1000];
  double SF_bkg_err[1000];
  double SF_edge[1000];
  double SF_edge_err[1000];
  double weight[1000];
  double AoE_ref = 0;
  double SF_dep_ref = 0;
  double SF_fep_ref = 0;
  double SF_sep_ref = 0;
  double SF_2614_ref = 0;
  double SF_bkg_ref = 0;
  double SF_edge_ref = 0;

  int n_lines = 0;
  while (1){
    file >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> x10 >> x11 >> x12 >> x13 >> x14 >> x15 >> x16;
    if (!file.good()) break;
    if ( x16 == 0 ){
      AoE_ref = x2;
      SF_dep_ref = x4;
      SF_fep_ref = x6;
      SF_sep_ref = x8;
      SF_2614_ref = x10;
      SF_bkg_ref = x12;
      SF_edge_ref = x14;
    }
    else {
      weight[n_lines] = x16;
      AoE[n_lines] = x2;
      AoE_err[n_lines] = x3;
      SF_dep[n_lines] = x4;
      SF_dep_err[n_lines] = x5;
      SF_fep[n_lines] = x6;
      SF_fep_err[n_lines] = x7;
      SF_sep[n_lines] = x8;
      SF_sep_err[n_lines] = x9;
      SF_2614[n_lines] = x10;
      SF_2614_err[n_lines] = x11;
      SF_bkg[n_lines] = x12;
      SF_bkg_err[n_lines] = x13;
      SF_edge[n_lines] = x14;
      SF_edge_err[n_lines] = x15;
      n_lines++;
    }
  }
  file.close(); 
  cout << "n. weight " << n_lines << endl;

  TCanvas *c1 = new TCanvas("c","cii",1000, 600);
  //c1->Divide(2,3);

  TLine *lin1 = new TLine(weight[0], AoE_ref,     weight[n_lines-1], AoE_ref);
  TLine *lin2 = new TLine(weight[0], SF_dep_ref,  weight[n_lines-1], SF_dep_ref);
  TLine *lin3 = new TLine(weight[0], SF_bkg_ref,  weight[n_lines-1], SF_bkg_ref);
  TLine *lin4 = new TLine(weight[0], SF_sep_ref,  weight[n_lines-1], SF_sep_ref);
  TLine *lin5 = new TLine(weight[0], SF_fep_ref,  weight[n_lines-1], SF_fep_ref);
  TLine *lin6 = new TLine(weight[0], SF_2614_ref, weight[n_lines-1], SF_2614_ref);
  TLine *lin7 = new TLine(weight[0], SF_edge_ref, weight[n_lines-1], SF_edge_ref);

  lin1->SetLineWidth(1);
  lin2->SetLineWidth(1);
  lin3->SetLineWidth(1);
  lin4->SetLineWidth(1);
  lin5->SetLineWidth(1);
  lin6->SetLineWidth(1);
  lin7->SetLineWidth(1);

  lin1->SetLineStyle(2);
  lin2->SetLineStyle(2);
  lin3->SetLineStyle(2);
  lin4->SetLineStyle(2);
  lin5->SetLineStyle(2);
  lin6->SetLineStyle(2);
  lin7->SetLineStyle(2);
 
  lin1->SetLineColor(4);
  lin2->SetLineColor(2);
  lin3->SetLineColor(kMagenta-7);
  lin4->SetLineColor(kMagenta+2);
  lin5->SetLineColor(8);
  lin6->SetLineColor(kTeal-7);
  lin7->SetLineColor(kCyan+1);

  TGraph *g_1 = new TGraphErrors(n_lines,weight, AoE, 0, AoE_err);
  
  TLegend *leg1 = new TLegend(0.4,0.8,0.6,0.9);
  leg1->AddEntry(g_1,"FWHM_{A/E}","PL");
  leg1->AddEntry(lin1,"Standard","l");

  g_1->SetTitle(Form("chn[%d]",chn));
  g_1->GetXaxis()->SetTitle("Weight");
  g_1->GetYaxis()->SetTitle("FWHM_{A/E} [%]");          

  g_1->SetMarkerColor(4);
  g_1->SetLineColor(4);
  g_1->SetLineWidth(2);
  g_1->SetMarkerStyle(20);

  //c1->cd(1);
  gPad->SetLogx();
  g_1->Draw("AP");
  leg1->Draw("same");
  lin1->Draw("same");
  c1->Update();
  c1->Print(Form("%s/FWHM_AoE_chn%d.pdf",resdir,chn));

  //TMultiGraph *mg = new TMultiGraph();
  TGraph *g_dep = new TGraphErrors(n_lines,weight, SF_dep, 0, SF_dep_err);
  TGraph *g_fep = new TGraphErrors(n_lines,weight, SF_fep, 0, SF_fep_err);
  TGraph *g_sep = new TGraphErrors(n_lines,weight, SF_sep, 0, SF_sep_err);
  TGraph *g_2614 = new TGraphErrors(n_lines,weight, SF_2614, 0,SF_2614_err);
  TGraph *g_bkg = new TGraphErrors(n_lines,weight, SF_bkg, 0,SF_bkg_err);
  TGraph *g_edge = new TGraphErrors(n_lines,weight, SF_edge, 0,SF_edge_err);  

  TLegend *leg2 = new TLegend(0.4,0.8,0.6,0.9);
  leg2->AddEntry(g_dep,"DEP","PL");
  leg2->AddEntry(lin2,"Standard","l");
  TLegend *leg3 = new TLegend(0.4,0.8,0.6,0.9);
  leg3->AddEntry(g_bkg,"ROI","PL");
  leg3->AddEntry(lin3,"Standard","l");
  TLegend *leg4 = new TLegend(0.4,0.8,0.6,0.9);
  leg4->AddEntry(g_sep,"SEP","PL");
  leg4->AddEntry(lin4,"Standard","l");
  TLegend *leg5 = new TLegend(0.4,0.8,0.6,0.9);
  leg5->AddEntry(g_fep,"FEP_{1621}","PL");
  leg5->AddEntry(lin5,"Standard","l");
  TLegend *leg6 = new TLegend(0.4,0.8,0.6,0.9);
  leg6->AddEntry(g_2614,"FEP_{2615}","PL");
  leg6->AddEntry(lin6,"Standard","l");
  TLegend *leg7 = new TLegend(0.4,0.8,0.6,0.9);
  leg7->AddEntry(g_edge,"C.Edge","PL");
  leg7->AddEntry(lin7,"Standard","l");

  g_dep->SetMarkerColor(2);
  g_fep->SetMarkerColor(8);
  g_sep->SetMarkerColor(kMagenta+2);
  g_2614->SetMarkerColor(kTeal-7);
  g_bkg->SetMarkerColor(kMagenta-7);
  g_edge->SetMarkerColor(kCyan+1);

  g_dep->SetLineColor(2);
  g_fep->SetLineColor(8);
  g_sep->SetLineColor(kMagenta+2);
  g_2614->SetLineColor(kTeal-7);
  g_bkg->SetLineColor(kMagenta-7);
  g_edge->SetLineColor(kCyan+1);

  g_dep->SetLineWidth(2);
  g_fep->SetLineWidth(2);
  g_sep->SetLineWidth(2);
  g_2614->SetLineWidth(2);
  g_bkg->SetLineWidth(2);
  g_edge->SetLineWidth(2);

  g_dep->SetMarkerStyle(20);
  g_fep->SetMarkerStyle(22);
  g_sep->SetMarkerStyle(23);
  g_2614->SetMarkerStyle(20);
  g_bkg->SetMarkerStyle(20);
  g_edge->SetMarkerStyle(20);

  //c1->cd(4);
  g_sep->SetTitle(Form("chn[%d]",chn));
  g_sep->GetXaxis()->SetTitle("Weight");
  g_sep->GetYaxis()->SetTitle("Survival Fraction [%]");
  g_sep->Draw("AP");
  gPad->SetLogx();
  leg4->Draw();
  lin4->Draw();
  c1->Update();
  c1->Print(Form("%s/SF_SEP_chn%d.pdf",resdir,chn));
  
  //c1->cd(5);
  g_fep->SetTitle(Form("chn[%d]",chn));
  g_fep->GetXaxis()->SetTitle("Weight");
  g_fep->GetYaxis()->SetTitle("Survival Fraction [%]");
  g_fep->Draw("AP");
  gPad->SetLogx();
  leg5->Draw("same");
  lin5->Draw();
  c1->Update();
  c1->Print(Form("%s/SF_FEP_chn%d.pdf",resdir,chn));
    
  //c1->cd(6);
  g_2614->SetTitle(Form("chn[%d]",chn));
  g_2614->GetXaxis()->SetTitle("Weight");
  g_2614->GetYaxis()->SetTitle("Survival Fraction [%]");
  g_2614->Draw("AP");
  gPad->SetLogx();
  leg6->Draw("same");
  lin6->Draw();
  c1->Update();
  c1->Print(Form("%s/SF_2614_chn%d.pdf",resdir,chn));
  
  //c1->cd(2);
  leg2->Draw("same");
  g_dep->SetTitle(Form("chn[%d]",chn));
  g_dep->GetXaxis()->SetTitle("Weight");
  g_dep->GetYaxis()->SetTitle("Survival Fraction [%]");
  g_dep->Draw("AP");
  leg2->Draw("same");
  gPad->SetLogx();
  lin2->Draw();
  c1->Update();
  c1->Print(Form("%s/SF_DEP_chn%d.pdf",resdir,chn));
  
  //c1->cd(3);
  g_bkg->GetXaxis()->SetTitle("Weight");
  g_bkg->GetYaxis()->SetTitle("Survival Fraction [%]");
  g_bkg->SetTitle(Form("chn[%d]",chn));
  g_bkg->Draw("AP");
  gPad->SetLogx();
  leg3->Draw("same");
  lin3->Draw();
  c1->Update();
  c1->Print(Form("%s/SF_bkg_chn%d.pdf",resdir,chn));
  //c1->Print(Form("%s/SF_FWHM_chn%d.pdf",resdir,chn));
  
  g_edge->GetXaxis()->SetTitle("Weight");
  g_edge->GetYaxis()->SetTitle("Survival Fraction [%]");
  g_edge->SetTitle(Form("chn[%d]",chn));
  g_edge->Draw("AP");
  gPad->SetLogx();
  leg7->Draw("same");
  lin7->Draw();
  c1->Update();
  c1->Print(Form("%s/SF_edge_chn%d.pdf",resdir,chn));

  // SAVE ON ROOT FILE
  TFile *out = new TFile(Form("%s/PlotSF_chn%d.root",resdir,chn),"RECREATE"); 
  out->cd();
  g_1->Write();
  g_dep->Write();
  g_fep->Write();
  g_sep->Write();
  g_bkg->Write();
  g_2614->Write();
  g_edge->Write();
  out->Close();
  
  //myapp->Run();
  return 0;
 }
