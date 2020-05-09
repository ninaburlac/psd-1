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
  char *text = argv[3];
  //memcpy(text, name + 12,5);
  
  ifstream file;
  file.open(name);
  
  double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
  double x11, x12, x13, x14;
  double chn[1000];
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
  
  int n_lines = 0;
  while (1){
    file >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> x10 >> x11 >> x12 >> x13 >> x14;
    if (!file.good()) break;
    chn[n_lines] = x1;
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
    n_lines++;
  }
  file.close(); 
  cout << "n. detector " << n_lines << endl;

  //--------------------- Draw the Survival Fraction-------------------------
  
  TCanvas *c1 = new TCanvas("c","SurvivalFraction");
  /*      
  TGraph *g_dep = new TGraphErrors(n_lines,chn, SF_dep, 0, SF_dep_err);
  TGraph *g_fep = new TGraphErrors(n_lines,chn, SF_fep, 0, SF_fep_err);
  TGraph *g_sep = new TGraphErrors(n_lines,chn, SF_sep, 0, SF_sep_err);
  TGraph *g_2614 = new TGraphErrors(n_lines,chn, SF_2614, 0,SF_2614_err);
  TGraph *g_bkg = new TGraphErrors(n_lines,chn, SF_bkg, 0,SF_bkg_err);
  */
  
  TGraph *g_dep = new TGraph(n_lines,chn, SF_dep);                            
  TGraph *g_fep = new TGraph(n_lines,chn, SF_fep);
  TGraph *g_sep = new TGraph(n_lines,chn, SF_sep);                            
  TGraph *g_2614 = new TGraph(n_lines,chn, SF_2614);                          
  TGraph *g_bkg = new TGraph(n_lines,chn, SF_bkg); 
  
  g_dep->SetTitle("run 98");
  g_dep->GetXaxis()->SetTitle("Detector number");
  g_dep->GetYaxis()->SetTitle("Survival Fraction [%]");
  g_dep->GetXaxis()->SetRangeUser(-0.5,40.5);
  g_dep->GetYaxis()->SetRangeUser(-0.5,100);

  g_dep->SetMarkerColor(2);
  g_fep->SetMarkerColor(8);
  g_sep->SetMarkerColor(kMagenta+2);
  g_2614->SetMarkerColor(kTeal-7);
  g_bkg->SetMarkerColor(kMagenta-7);

  g_dep->SetLineColor(2);
  g_fep->SetLineColor(8);
  g_sep->SetLineColor(kMagenta+2);
  g_2614->SetLineColor(kTeal-7);
  g_bkg->SetLineColor(kMagenta-7);

  g_dep->SetLineWidth(2);
  g_fep->SetLineWidth(2);
  g_sep->SetLineWidth(2);
  g_2614->SetLineWidth(2);
  g_bkg->SetLineWidth(2);

  g_dep->SetMarkerStyle(20);
  g_fep->SetMarkerStyle(22);
  g_sep->SetMarkerStyle(23);
  g_2614->SetMarkerStyle(20);
  g_bkg->SetMarkerStyle(20);

  g_dep->Draw("AP");
  g_fep->Draw("P");
  g_sep->Draw("P");
  g_2614->Draw("P");
  g_bkg->Draw("P");
  c1->Update();
  
  TText *xlab = new TText();
  xlab-> SetNDC();
  xlab -> SetTextFont(1);
  xlab -> SetTextColor(1);
  xlab -> SetTextSize(0.03);
  xlab -> SetTextAlign(22);
  xlab -> SetTextAngle(0);
  xlab -> DrawText(0.57, 0.92, "String 1                                     String 3                      String 4                                    String 6           ");
  c1->Update();

  TLine *lin1 = new TLine(7.5,c1->GetUymin(),7.5,c1->GetUymax());
  TLine *lin2 = new TLine(10.5,c1->GetUymin(),10.5,c1->GetUymax());
  TLine *lin3 = new TLine(18.5,c1->GetUymin(),18.5,c1->GetUymax());
  TLine *lin4 = new TLine(26.5,c1->GetUymin(),26.5,c1->GetUymax());
  TLine *lin5 = new TLine(29.5,c1->GetUymin(),29.5,c1->GetUymax());
  TLine *lin16 = new TLine(35.5,c1->GetUymin(),35.5,c1->GetUymax());
  
  lin1->SetLineWidth(1);
  lin2->SetLineWidth(1);
  lin3->SetLineWidth(1); 
  lin4->SetLineWidth(1);
  lin5->SetLineWidth(1);
  lin16->SetLineWidth(1);

  lin1->SetLineStyle(2);
  lin2->SetLineStyle(2);
  lin3->SetLineStyle(2);
  lin4->SetLineStyle(2);
  lin5->SetLineStyle(2);
  lin16->SetLineStyle(2);

  lin1->Draw();
  lin2->Draw();
  lin3->Draw(); 
  lin4->Draw();
  lin5->Draw();
  lin16->Draw();
  c1->Update();

  TLegend *leg = new TLegend(0.82,0.58,0.95,0.77);
  leg->AddEntry(g_dep,"DEP","PL");
  leg->AddEntry(g_fep,"FEP_{1621}","PL");
  leg->AddEntry(g_sep,"SEP","PL");
  leg->AddEntry(g_2614,"FEP_{2615}","PL");
  leg->AddEntry(g_bkg,"ROI","PL");
  leg->Draw();
  c1->Update();
  c1->Print(Form("%s/SF_peaks_%s.pdf",resdir,text));
  
  ////////////////////////////////////////////////////////////////
  //---------------- Draw the FWHM of AoE------------------------
  ///////////////////////////////////////////////////////////////

  TCanvas *c2 = new TCanvas("c2","AoE");
  //TGraph *g_1 = new TGraphErrors(n_lines,chn,AoE,0,AoE_err);
  TGraph *g_1 = new TGraph(n_lines,chn,AoE);
  
  g_1->SetTitle("run 98");
  g_1->GetXaxis()->SetTitle("Detector number");
  g_1->GetYaxis()->SetTitle("FWHM_{A/E} [%]");
  g_1->GetXaxis()->SetTitleSize(0.04);
  g_1->GetYaxis()->SetTitleSize(0.04);
  g_1->GetXaxis()->SetRangeUser(-0.5, 40.5);
  //g_1->GetYaxis()->SetRangeUser(1.3, 4.5);
  
  g_1->SetMarkerColor(2);
  g_1->SetLineColor(2);
  g_1->SetLineWidth(2);
  g_1->SetMarkerStyle(20);
  
  g_1->Draw("AP");
  c2->Update();

  TText *xlabel = new TText();
  xlabel-> SetNDC();
  xlabel -> SetTextFont(1);
  xlabel -> SetTextColor(1);
  xlabel -> SetTextSize(0.03);
  xlabel -> SetTextAlign(22);
  xlabel -> SetTextAngle(0);
  xlabel -> DrawText(0.57, 0.92, "String 1                                     String 3                      String 4                                    String 6           ");
  c2->Update();

  TLine *line1 = new TLine(7.5,c2->GetUymin(),7.5,c2->GetUymax());
  TLine *line2 = new TLine(10.5,c2->GetUymin(),10.5,c2->GetUymax());
  TLine *line3 = new TLine(18.5,c2->GetUymin(),18.5,c2->GetUymax());                               
  TLine *line4 = new TLine(26.5,c2->GetUymin(),26.5,c2->GetUymax());
  TLine *line5 = new TLine(29.5,c2->GetUymin(),29.5,c2->GetUymax());
  TLine *line6 = new TLine(35.5,c2->GetUymin(),35.5,c2->GetUymax());  
  
  line1->SetLineWidth(1);
  line2->SetLineWidth(1);
  line3->SetLineWidth(1);
  line4->SetLineWidth(1);
  line5->SetLineWidth(1);
  line6->SetLineWidth(1); 

  line1->SetLineStyle(2);
  line2->SetLineStyle(2);
  line3->SetLineStyle(2);                                                         
  line4->SetLineStyle(2);
  line5->SetLineStyle(2);
  line6->SetLineStyle(2);  

  line1->Draw();
  line2->Draw();
  line3->Draw();        
  line4->Draw();
  line5->Draw();
  line6->Draw();
  c2->Update();
  c2->Print(Form("%s/FWHM_AoE_%s.pdf",resdir,text));
  
  // SAVE ON ROOT FILE
  TFile *out = new TFile(Form("%s/PlotSF.root",resdir),"RECREATE"); 
  out->cd();
  c1->Write();
  c2->Write();
  out->Close();
  
  //myapp->Run();
  return 0;
}
